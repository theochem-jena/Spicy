{-|
Module      : Spicy.Parser
Description : Parsers for chemical data formats and computational chemistry output files.
Copyright   : Phillip Seeber, 2019
License     : GPL-3
Maintainer  : phillip.seeber@uni-jena.de
Stability   : experimental
Portability : POSIX, Windows

-}
{-# LANGUAGE OverloadedStrings #-}
module Spicy.Parser
( -- * Chemical Data Formats
  parseXYZ
, parseTXYZ
, parseMOL2
, parsePDB
, parseSpicy
-- * Generic formats
, parseHMatrix
) where
import           Control.Applicative
import qualified Data.Array.Accelerate     as A
import           Data.Attoparsec.Text.Lazy
import           Data.Char
import           Data.Either
import           Data.Foldable
import           Data.IntMap               (IntMap)
import qualified Data.IntMap               as IM
import           Data.IntSet               (IntSet)
import qualified Data.IntSet               as IS
import           Data.Maybe
import           Data.Sequence             (Seq (..))
import qualified Data.Sequence             as S
import qualified Data.Text                 as TS
import qualified Data.Text.Lazy            as TL
import qualified Data.Text.Lazy.Read       as TL
import           Data.Tuple
import           Lens.Micro.Platform
import           Prelude                   hiding (cycle, foldl1, foldr1, head,
                                            init, last, maximum, minimum, tail,
                                            take, takeWhile, (!!))
import           Spicy.Molecule.Util
import           Spicy.Types
import           Text.Read
-- import Data.Monoid
-- import Data.Sequence (Seq)
-- import Debug.Trace
-- import Control.DeepSeq


{-|
Make a parser optional and wrap it in a 'Maybe'.
-}
maybeOption :: Parser a -> Parser (Maybe a)
maybeOption p = option Nothing (Just <$> p)

{-|
Parser for skiping non line breaking space.
-}
skipSpace' :: Parser ()
skipSpace' = do
  _ <- takeWhile (`elem` [' ', '\t', '\f', '\v'])
  return ()

{-|
Convert strict text to lazy text.
-}
textS2L :: TS.Text -> TL.Text
textS2L = TL.pack . TS.unpack

----------------------------------------------------------------------------------------------------
{-|
Parse a .xyz file (has no connectivity, atom types or partioal charges).
-}
parseXYZ :: Parser Molecule
parseXYZ = do
  nAtoms <- skipSpace' *> decimal
  label  <- skipSpace *> takeWhile (not . isEndOfLine) <* skipSpace
  atoms  <- count nAtoms xyzLineParser
  return Molecule
    { _molecule_Label    = textS2L label
    , _molecule_Atoms    = IM.fromList $ zip [ 0 .. ] atoms
    , _molecule_Bonds    = IM.empty
    , _molecule_SubMol   = S.empty
    , _molecule_Energy   = Nothing
    , _molecule_Gradient = Nothing
    , _molecule_Hessian  = Nothing
    }
  where
    xyzLineParser :: Parser Atom
    xyzLineParser = do
      cElement <- skipSpace' *> many1 letter
      x        <- skipSpace' *> double
      y        <- skipSpace' *> double
      z        <- skipSpace' *> double
      skipSpace
      return Atom
        { _atom_Element     = fromMaybe H . readMaybe $ cElement
        , _atom_Label       = ""
        , _atom_IsPseudo    = False
        , _atom_FFType      = ""
        , _atom_PCharge     = Nothing
        , _atom_Coordinates = S.fromList [x, y, z]
        }

{-|
Parse a Tinker XYZ formatted file. It has coordinates and might have connectivity and atom types.
This format and therefore parser are not using any layers (recursions of 'Molecule').
-}
parseTXYZ :: Parser Molecule
parseTXYZ = do
  _nAtoms     <- skipSpace' *> (decimal :: Parser Int)
  label       <- skipSpace' *> takeWhile (not . isEndOfLine) <* skipSpace
  conAndAtoms <- many1 txyzLineParser
  return Molecule
    { _molecule_Label    = textS2L label
    , _molecule_Atoms    = IM.fromList . map (\a -> (a ^. _1 . _1, a ^. _2)) $ conAndAtoms
    , _molecule_Bonds    = IM.fromList . map fst $ conAndAtoms
    , _molecule_SubMol   = S.empty
    , _molecule_Energy   = Nothing
    , _molecule_Gradient = Nothing
    , _molecule_Hessian  = Nothing
    }
  where
    -- Parsing a single line of atoms. Tinker's format keeps bonds associated with atoms. So a tuple
    -- suitable to construct the 'IntMap' is returned additional to the pure atoms.
    txyzLineParser :: Parser ((Int, IntSet), Atom)
    txyzLineParser = do
      index           <- skipSpace' *> ((\a -> a - 1) <$> decimal)
      cElement        <- skipSpace' *> many1 letter
      x               <- skipSpace' *> double
      y               <- skipSpace' *> double
      z               <- skipSpace' *> double
      mFFType         <- skipSpace' *> maybeOption decimal
      connectivityRaw <- skipSpace' *> many' columnDecimal
      endOfLine
      return
        ( (index, IS.fromList connectivityRaw)
        , Atom
            { _atom_Element      = fromMaybe H . readMaybe $ cElement
            , _atom_Label        = ""
            , _atom_IsPseudo     = False
            , _atom_FFType       =
                case mFFType of
                  Nothing -> ""
                  Just a  -> TL.pack . show $ (a :: Int)
            , _atom_PCharge      = Nothing
            , _atom_Coordinates  = S.fromList [x, y, z]
            }
          )
    -- Parse multiple non-line-breaking whitespace separated decimals.
    columnDecimal :: Parser Int
    columnDecimal = skipSpace' *> decimal <* skipSpace'

{-|
Parse the "interesting" fields of a MOL2 file. This contains partial charges as well as
connectivity. There is no special understanding for the atom types, that are available in MOL2
files. They will simply be treated as the force field string. See
<http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf>.
-}
parseMOL2 :: Parser Molecule
parseMOL2 = do
  (label, nAtoms, nBonds) <- moleculeParser
  atomsLabeled            <- atomParser nAtoms
  bonds                   <- bondParser nBonds
  let -- Construct the top layer of the molecule.
      atoms        = IM.fromList . toList . fmap (\(ind, _, atom) -> (ind, atom)) $ atomsLabeled
      -- From the groups of same substructure ID fragments, build molecules.
      subMols = makeSubMolsFromAnnoAtoms atomsLabeled bonds
  return Molecule
    { _molecule_Label    = label
    , _molecule_Atoms    = atoms
    , _molecule_Bonds    = bonds
    , _molecule_SubMol   = subMols
    , _molecule_Energy   = Nothing
    , _molecule_Gradient = Nothing
    , _molecule_Hessian  = Nothing
    }
  where
    -- Parse the @<TRIPOS>MOLECULE block of MOL2.
    moleculeParser :: Parser (TL.Text, Int, Maybe Int)
    moleculeParser = do
      _header     <- manyTill anyChar (string "@<TRIPOS>MOLECULE") <* endOfLine
      -- Line 1 -> "mol_name"
      label       <- takeWhile (not . isEndOfLine) <* endOfLine
      -- Line 2 -> "num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]"
      nAtoms      <- skipSpace' *> decimal
      nBonds      <- maybeOption $ skipSpace' *> decimal
      _nSubMols   <- maybeOption $ skipSpace' *> (decimal :: Parser Int)
      _nFeatures  <- maybeOption $ skipSpace' *> (decimal :: Parser Int)
      _nSets      <- maybeOption $ skipSpace' *> (decimal :: Parser Int) <* skipSpace' <* endOfLine
      -- Line 3 -> "mol_type"
      _molType    <-
            skipSpace'
        *>  string "SMALL"
        <|> string "BIOPOLYMER"
        <|> string "PROTEIN"
        <|> string "NUCLEIC_ACID"
        <|> string "SACCHARIDE"
        <*  skipSpace
      -- Line 4 -> "charge_type"
      _chargeType <-
            skipSpace
        *>( string "NO_CHARGES"
        <|> string "DEL_RE"
        <|> string "GASTEIGER"
        <|> string "GAST_HUCK"
        <|> string "HUCKEL"
        <|> string "PULLMAN"
        <|> string "GAUSS80_CHARGES"
        <|> string "AMPAC_CHARGES"
        <|> string "MULLIKEN_CHARGES"
        <|> string "DICT_CHARGES"
        <|> string "MMFF94_CHARGES"
        <|> string "USER_CHARGES"
        )<* skipSpace'
        <*  endOfLine
      -- Line 5 -> "[status_bits"
      _statusBit  <- maybeOption $ skipSpace' *> many1 letter <* skipSpace
      -- Line 6 -> "[mol_comment]]"
      _comment    <- maybeOption $ skipSpace' *> takeWhile (not . isEndOfLine) <* skipSpace
      return (textS2L label, nAtoms, nBonds)
    --
    -- Parse the @<TRIPOS>ATOM block of MOL2. This will give:
    -- (index of the atom, (substructure ID, substructure name), atom).
    atomParser :: Int -> Parser (Seq (Int, (Int, TL.Text), Atom))
    atomParser nAtoms = do
      _header <- manyTill anyChar (string "@<TRIPOS>ATOM") <* endOfLine
      -- Parse multiple lines of ATOM data.
      atoms   <- count nAtoms $ do -- atomLineParser
        index    <- skipSpace' *> ((\a -> a - 1) <$> decimal)
        -- Most often this will be the element symbol.
        label    <- skipSpace' *> takeWhile (not . isHorizontalSpace)
        -- x, y and z coordinates
        x        <- skipSpace' *> double
        y        <- skipSpace' *> double
        z        <- skipSpace' *> double
        -- Parse the chemical element, which is actually the first part of the SYBYL atom type.
        cElem    <- skipSpace' *> many1 letter
        -- A dot often separates the element from the type of this element.
        ffdot    <- maybeOption $ char '.'
        -- And after the dot the rest of the SYBYL atom type might come.
        ffType   <- maybeOption $ takeWhile (not . isHorizontalSpace)
        -- The substructure ID. This is the identifier to identify sub molecules.
        subID    <- skipSpace' *> (decimal :: Parser Int)
        -- The substructure Name. This should be used as label for the sub molecule.
        subName  <- skipSpace' *> takeWhile (not . isHorizontalSpace)
        -- The partial charge
        pCharge  <- skipSpace' *> double <* skipSpace
        return
          ( index
          , (subID, textS2L subName)
          , Atom
              { _atom_Element     = fromMaybe H . readMaybe $ cElem
              , _atom_Label       = textS2L label
              , _atom_IsPseudo    = False
              , _atom_FFType      =
                  (TL.pack cElem)
                  `TL.append`
                  ( TL.pack . (\c -> case c of
                      Just a  -> [a]
                      Nothing -> ""
                    ) $ ffdot
                  )
                  `TL.append`
                  (fromMaybe "" $ textS2L <$> ffType)
              , _atom_PCharge     = Just pCharge
              , _atom_Coordinates = S.fromList [x, y, z]
              }
          )
      return $ S.fromList atoms
    --
    -- Parse the @<TRIPOS>BOND part. Unfortunately, the bonds in the MOL2 format are unidirectiorial
    -- and need to be flipped to.
    bondParser :: Maybe Int -> Parser (IntMap IntSet)
    bondParser nBonds = do
      let -- How often to parse bond fields depends on if the number of bonds has been specified.
          nParser =
            case nBonds of
              Nothing -> many'
              Just n  -> count n
      _header  <- manyTill anyChar (string "@<TRIPOS>BOND") <* endOfLine
      uniBonds <- nParser $ do
        -- Bond id, which does not matter.
        _id    <- skipSpace' *> (decimal :: Parser Int)
        -- Origin atom index
        origin <- skipSpace' *> ((\a -> a - 1) <$> decimal :: Parser Int)
        -- Target atom index
        target <- skipSpace' *> ((\a -> a - 1) <$> decimal :: Parser Int)
        -- Bond type, which we don't care about.
        _type  <- skipSpace' *> takeWhile (not . isSpace) <* skipSpace
        return (origin, target)
      let -- Make the bonds bidirectorial
          bondTupleSeq = S.fromList uniBonds
          bondsForth   = groupTupleSeq bondTupleSeq
          bondsBack    = groupTupleSeq $ swap <$> bondTupleSeq
          bonds        = bondsForth <> bondsBack
      return bonds

{-|
Parse a PDB file as described in
<ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf>.
-}
parsePDB :: Parser Molecule
parsePDB = do
  -- Parse atoms only and ignore other fiels
  atomsLabeled <- S.fromList <$> many1 atomParser
  bonds        <- IM.fromList <$> many' connectParser
  -- links <- undefined --many' linkParser
  let -- Transform the informations from the parsers.
      atomsIM = IM.fromList . toList .  fmap (\(ind, _, atom) -> (ind, atom)) $ atomsLabeled
      subMols = makeSubMolsFromAnnoAtoms atomsLabeled bonds
  return Molecule
    { _molecule_Label    = ""
    , _molecule_Atoms    = atomsIM
    , _molecule_Bonds    = bonds
    , _molecule_SubMol   = subMols
    , _molecule_Energy   = Nothing
    , _molecule_Gradient = Nothing
    , _molecule_Hessian  = Nothing
    }
  where
    -- This parser works a little bit different than the others, as this is a fixed columnd witdth
    -- format and we can't rely on white spaces or fields really containing a value. A tuple of
    -- following structure is returned: (Index of atom, (subMolID, subMolName), atom)
    atomParser :: Parser (Int, (Int, TL.Text), Atom)
    atomParser = do
      -- First check, that this line really starts with an ATOM or HETATM record (6 characters).
      -- Columns 1-6: record type.
      _recordStart <- TL.pack <$> manyTill anyChar (string "ATOM  " <|> string "HETATM")
      -- Then take the rest of the line, till the end of line is reached.
      recordRest   <- textS2L <$> takeWhile (not . isEndOfLine) <* endOfLine
      let -- Recombine the line and split it according to the PDB format specifier.
          atomLine          = "ATOM  " `TL.append` recordRest
          -- Now according to the PDB specification. Fields exaclty named as in the PDF with "c" as
          -- prefix for chunk.
          -- Columns 1-6: record type.
          (_cAtom, rest1)       = TL.splitAt 6 atomLine
          -- Column 7-11: atom serial number
          (cSerial, rest2)      = TL.splitAt 5 rest1
          -- Column 13-16: atom name
          (cName, rest3)        = TL.splitAt 4 . TL.drop 1 $ rest2
          -- Column 17: alternate location indicator
          (_cAltLoc, rest4)     = TL.splitAt 1 rest3
          -- Columnt 18-20: residue name (use this as submolecule label)
          (cResName, rest5)     = TL.splitAt 3 . TL.drop 1 $ rest4
          -- Column 22: chain identifier
          (_cChainID, rest6)    = TL.splitAt 1 rest5
          -- Column 23-26: residue sequence number
          (cResSeq, rest7)      = TL.splitAt 4 rest6
          -- Column 27: Code for insertion of residue
          (_cICode, rest8)      = TL.splitAt 1 rest7
          -- Column 31-38: Orthogonal coordinates for x in Angstrom
          (cX, rest9)           = TL.splitAt 8 . TL.drop 3 $ rest8
          -- Column 39-46: Orthogonal coordinates for y in Angstrom
          (cY, rest10)          = TL.splitAt 8 rest9
          -- Column 47-54: Orthogonal coordinates for z in Angstrom
          (cZ, rest11)          = TL.splitAt 8 rest10
          -- Column 55-60: Occupancy
          (_cOccupancy, rest12) = TL.splitAt 6 rest11
          -- Column 61-66: temperature factor
          (_cTempFactor, rest13)= TL.splitAt 6 . TL.drop 10 $ rest12
          -- Column 77-78: element
          (cElement, rest14)    = TL.splitAt 2 rest13
          -- Columnt 79-80: charge on the atom
          (cCharge, _rest15)    = TL.splitAt 2 rest14
          --
          -- Now parse all fields of interest using text readers.
          aSerial      = fst <$> (TL.decimal . TL.strip $ cSerial)
          aResSeq      = fst <$> (TL.decimal . TL.strip $ cResSeq)
          aElement     =
            let elemMaybe = (readMaybe :: String -> Maybe Element) . TL.unpack . TL.strip $ cElement
            in  case elemMaybe of
                  Nothing -> Left "parsePDB: Could not read the element symbol."
                  Just e  -> Right e
          aLabel       = TL.strip cName
          aFFType      = ""
          aPCharge     =
            let pChargeMaybe = fst <$> (TL.double . TL.strip $ cCharge)
            in  case pChargeMaybe of
                  Left _   -> Nothing
                  Right pC -> Just pC
          aCoordinates = traverse (fmap fst . TL.double . TL.strip) . S.fromList $ [cX, cY, cZ]
          -- Use the Atom constructor applicatively to get Either String Atom. This relies on the
          -- order of the atom arguments.
          eitherAtom   =
            Atom
            <$> aElement            -- _atom_Element
            <*> pure aLabel         -- _atom_Label
            <*> pure False          -- _atom_IsPseudo
            <*> pure aFFType        -- _atom_FFType
            <*> pure aPCharge       -- _atom_FFType
            <*> aCoordinates        -- _atom_Coordinates
      -- If all the non-Attoparsec parse actions succeeded, return an Attoparsec result or fail with
      -- an Attoparsec error.
      case (aSerial, (aResSeq, cResName), eitherAtom) of
        (Right ind, (Right sInd, _), Right a) -> return (ind, (sInd, cResName), a)
        (Left indErr, _, _)                   -> fail indErr
        (_, (Left sIndErr, _), _)             -> fail sIndErr
        (_, _, Left aErr)                     -> fail aErr
    --
    -- Parse CONECT fields of the PDB. PDB bonds are bidirectorial, so no swapping required.
    connectParser :: Parser (Int, IntSet)
    connectParser = do
      -- First check, that this line really starts with a CONECT record (6 characters).
      -- Columns 1-6: record type.
      _recordStart <- TL.pack <$> manyTill anyChar (string "CONECT")
      -- Then take the rest of the line, till the end of line is reached.
      recordRest   <- textS2L <$> takeWhile (not . isEndOfLine) <* endOfLine
      let -- Recombine the line and split them according to PDB standard.
          conectLine = "CONECT" `TL.append` recordRest
          -- Columns 1-6: "CONECT"
          (_cConect, rest1)  = TL.splitAt 6 conectLine
          -- Columns 7-11: serial of origin atom.
          (cOrigin, rest2)   = TL.splitAt 5 rest1
          -- Columns 12-16, 17-21, 22-26,27-31: serial of a target atoms.
          (cTarget1, rest3)  = TL.splitAt 5 rest2
          (cTarget2, rest4)  = TL.splitAt 5 rest3
          (cTarget3, rest5)  = TL.splitAt 5 rest4
          (cTarget4, _rest6) = TL.splitAt 5 rest5
          --
          -- Now parse all fields of interest using text readers.
          pOrigin  = fst <$> (TL.decimal . TL.strip $ cOrigin)
          pTargets :: Either String IntSet
          pTargets =
              fmap IS.fromList
            . sequence
            . filter isRight
            . map (fmap fst . TL.decimal . TL.strip)
            $ [cTarget1, cTarget2, cTarget3, cTarget4]
      case (pOrigin, pTargets) of
        (Right o, Right t) -> return (o, t)
        (Left oErr, _)     -> fail oErr
        (_, Left tErr)     -> fail tErr

{-|
Parser for the Spicy format used in this program. Represents fully all informations stored in the
'Molecule' type.
-}
parseSpicy :: Parser Molecule
parseSpicy = undefined {- do
  _ <- string "#Spicy-Format v0.2"
  endOfLine
  skipSpace
  _ <- string "#Spicy-Molecule"
  endOfLine
  _ <- string "  Label:"
  endOfLine
  _ <- string "    "
  mLabel <- takeTill isEndOfLine
  skipSpace
  mEnergy <- maybeOption parseEnergy
  mGradient <- maybeOption parseGradient
  mHessian <- maybeOption parseHessian
  skipSpace
  _ <- string "#Spicy-Atoms"
  skipSpace
  mAtoms <- many1 parseAtoms
  return Molecule
    { _molecule_Label    = TS.unpack mLabel
    , _molecule_Atoms    = R.fromList (R.Z R.:. length mAtoms) mAtoms
    , _molecule_Energy   = mEnergy
    , _molecule_Gradient = mGradient
    , _molecule_Hessian  = mHessian
    }
  where
    parseEnergy = do
      _ <- manyTill anyChar (string "Energy / Hartree:")
      skipSpace
      energy <- double
      return energy
    parseGradient = do
      _ <- manyTill anyChar (string "Gradient / Hartee/Bohr:")
      skipSpace
      gradient <- many1 $ do
        skipSpace
        gVal <- double
        return gVal
      return $ R.fromListUnboxed (R.Z R.:. (length gradient)) gradient
    parseHessian = do
      _ <- manyTill anyChar (string "Hessian / a.u.:")
      skipSpace
      hessian <- parseHMatrix
      return hessian
    parseAtoms = do
      -- indendation
      _ <- many' (char ' ' <|> char '\t')
      -- chemical element
      cElement <- many1 letter
      --
      _ <- count 4 (char ' ')
      -- label of the atom
      label <- count 6 anyChar

      --
      _ <- count 4 (char ' ')
      -- pseudo label
      pseudo <- anyChar
      --
      _ <- count 4 (char ' ')
      -- FFType
      ffType <- count 6 anyChar
      pChargeTest <- maybeOption $ do
        skipSpace
        (string "No")
      pCharge <- if pChargeTest == Nothing
        then Just <$> do
          _ <- many' (char ' ' <|> char '\t')
          double
        else return Nothing
      _ <- many' (char ' ' <|> char '\t')
      coordVec <- count 3 $ do
        coordComponent <- double
        _ <- many' (char ' ' <|> char '\t')
        return coordComponent
      _ <- many' (char ' ' <|> char '\t')
      connectivity <- many' $ do
        conAtom <- decimal
        _ <- many' (char ' ' <|> char '\t')
        return conAtom
      endOfLine
      return Atom
        { _atom_Element      = read cElement
        , _atom_Label        = T.unpack . T.strip . T.pack $ label
        , _atom_IsPseudo     = if pseudo == 'P' then True else False
        , _atom_FFType       = T.unpack . T.strip . T.pack $ ffType
        , _atom_PCharge      = pCharge
        , _atom_Coordinates  = R.fromListUnboxed (R.Z R.:. 3) coordVec
        , _atom_Connectivity = I.fromList connectivity
        }
-}

----------------------------------------------------------------------------------------------------
{-|
Parse the "show" instance output for HMatrix' matrix Type, including the dimension infos
-}
parseHMatrix :: Parser (A.Matrix Double)
parseHMatrix = undefined {- do
  skipSpace
  _ <- char '('
  dimX <- decimal
  _ <- char '>'
  _ <- char '<'
  dimY <- decimal
  _ <- char ')'
  skipSpace
  _ <- char '['
  _ <- many' (char ' ' <|> char '\t')
  rows <- do
    count dimX $ do
      singleRow <- count dimY parseElement
      skipSpace
      return singleRow
  _ <- char ']'
  --return $ fromRows . map fromList $ rows
  return $ R.fromListUnboxed (R.Z R.:. dimX R.:. dimY) . concat $ rows
  where
    parseElement :: Parser Double
    parseElement = do
      _ <- many' (char ' ' <|> char '\t')
      _ <- option ',' (char ',')
      _ <- many' (char ' ' <|> char '\t')
      element <- double
      _ <- many' (char ' ' <|> char '\t')
      _ <- option ',' (char ',')
      _ <- many' (char ' ' <|> char '\t')
      return element
-}
