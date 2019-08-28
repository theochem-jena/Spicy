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
import           Data.Tuple
import           Lens.Micro.Platform
import           Prelude                   hiding (cycle, foldl1, foldr1, head,
                                            init, last, maximum, minimum, tail,
                                            take, takeWhile, (!!))
import           Spicy.Molecule.Util
import           Spicy.Types
import           Text.Read
--import Data.Monoid
--import Data.Sequence (Seq)


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

----------------------------------------------------------------------------------------------------
{-|
Parse a .xyz file (has no connectivity, atom types or partioal charges).
-}
parseXYZ :: Parser Molecule
parseXYZ = do
  nAtoms <- skipSpace' *> decimal
  label  <- skipSpace *> takeWhile (/= '\n') <* skipSpace
  atoms  <- count nAtoms xyzLineParser
  return Molecule
    { _molecule_Label    = TL.pack . TS.unpack $ label
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
  label       <- skipSpace' *> takeWhile (/= '\n') <* skipSpace
  conAndAtoms <- many1 txyzLineParser
  return Molecule
    { _molecule_Label    = TL.pack . TS.unpack $ label
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
      atoms        = IM.fromList . map (\(ind, _, atom) -> (ind, atom)) $ atomsLabeled
      -- Sort the labeled atoms based on their substructure IDs.
      subMolGroups :: Seq (Seq (Int, (Int, TL.Text), Atom))
      subMolGroups =
          groupBy (\x y -> x ^. _2 . _1 == y ^. _2 . _1)
        . S.unstableSortOn (\(_, (subID, _), _) -> subID)
        . S.fromList
        $ atomsLabeled
      -- From the groups of same substructure ID fragments, build molecules.
      subMols = fmap (\g -> makeSubMolFromGroup g bonds) subMolGroups
      -- Assume we get properly sorted groups (all subIDs are them same): Build a new molecule from
      -- this group.
      makeSubMolFromGroup :: Seq (Int, (Int, TL.Text), Atom) -> IntMap IntSet -> Molecule
      makeSubMolFromGroup group bonds' =
        let atoms'      = IM.fromList . toList . fmap (\(ind, _, atom) -> (ind, atom)) $ group
            label'      = fromMaybe "" . (S.!? 0) . fmap (^. _2 . _2) $ group
            atomInds    = IM.keysSet atoms'
            -- Remove all bonds from the IntMap, that have origin on atoms not in this set and all
            -- target atoms that are not in the set.
            bondsCleaned :: IntMap IntSet
            bondsCleaned =
                IM.map (IS.filter (`IS.member` atomInds))
              $ bonds' `IM.restrictKeys` atomInds
        in  Molecule
              { _molecule_Label    = label'
              , _molecule_Atoms    = atoms'
              , _molecule_Bonds    = bondsCleaned
              , _molecule_SubMol   = S.empty
              , _molecule_Energy   = Nothing
              , _molecule_Gradient = Nothing
              , _molecule_Hessian  = Nothing
              }
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
      label       <- takeWhile (/= '\n') <* endOfLine
      --traceShowM label
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
      _comment    <- maybeOption $ skipSpace' *> takeWhile (/= '\n') <* skipSpace
      return (TL.pack . TS.unpack $ label, nAtoms, nBonds)
    --
    -- Parse the @<TRIPOS>ATOM block of MOL2. This will give:
    -- (index of the atom, (substructure ID, substructure name), atom).
    atomParser :: Int -> Parser [(Int, (Int, TL.Text), Atom)]
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
          , (subID, TL.pack . TS.unpack $ subName)
          , Atom
              { _atom_Element     = fromMaybe H . readMaybe $ cElem
              , _atom_Label       = TL.pack . TS.unpack $ label
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
                  (fromMaybe "" $ TL.pack . TS.unpack <$> ffType)
              , _atom_PCharge     = Just pCharge
              , _atom_Coordinates = S.fromList [x, y, z]
              }
          )
      return atoms
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

parsePDB :: Parser Molecule
parsePDB = undefined

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
