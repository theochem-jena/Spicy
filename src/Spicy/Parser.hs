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
import qualified Data.Array.Accelerate     as A
import           Data.Attoparsec.Text.Lazy
import qualified Data.IntMap               as IM
import           Data.Maybe
import qualified Data.Sequence             as S
import qualified Data.Text                 as TS
import qualified Data.Text.Lazy            as TL
import           Prelude                   hiding (cycle, foldl1, foldr1, head,
                                            init, last, maximum, minimum, tail,
                                            take, takeWhile, (!!))
import           Spicy.Types
import           Text.Read


{-|
Make a parser optional and wrap it in a 'Maybe'.
-}
maybeOption :: Parser a -> Parser (Maybe a)
maybeOption p = option Nothing (Just <$> p)

----------------------------------------------------------------------------------------------------
{-|
Parse a .xyz file (has no connectivity, atom types or partioal charges).
-}
parseXYZ :: Parser Molecule
parseXYZ = do
  nAtoms <- skipSpace *> decimal
  label  <- skipSpace *> takeWhile (/= '\n')
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
      cElement <- skipSpace *> many1 letter
      x        <- skipSpace *> double
      y        <- skipSpace *> double
      z        <- skipSpace *> double
      skipSpace
      return Atom
        { _atom_Element      = fromMaybe H . readMaybe $ cElement
        , _atom_Label        = ""
        , _atom_IsPseudo     = False
        , _atom_FFType       = ""
        , _atom_PCharge      = Nothing
        , _atom_Coordinates  = S.fromList [x, y, z]
        }


{-|
Parse a Tinker XYZ formatted file. It has coordinates and might have connectivity and atom types.
This format and therefore parser are not using any layers (recursions of 'Molecule').
-}
parseTXYZ :: Parser Molecule
parseTXYZ = undefined {-do
  skipSpace
  nAtoms <- decimal
  _ <- many' (char ' ' <|> char '\t')
  comment <- manyTill anyChar endOfLine
  conAndAtoms <- many1 txyzLineParser
  return Molecule
    { _molecule_Label    = comment
    , _molecule_Atoms    = VB.fromList . map snd $ conAndAtoms
    , _molecule_Bonds    = IM.fromList . map fst $ conAndAtoms
    , _molecule_Energy   = Nothing
    , _molecule_Gradient = Nothing
    , _molecule_Hessian  = Nothing
    }
  where
    -- Parsing a single line of atoms. Tinker's format keeps bonds associated with atoms. So a tuple
    -- suitable to construct the 'IntMap' is returned additional to the pure atoms.
    txyzLineParser :: Parser ((Int, IntSet), Atom)
    txyzLineParser = do
      skipSpace
      index <- (\a -> a - 1) <$> decimal
      skipSpace
      cElement <- many1 letter
      skipSpace
      x <- double
      skipSpace
      y <- double
      skipSpace
      z <- double
      skipSpace
      mFFType <- maybeOption decimal
      _ <- many' (char ' ' <|> char '\t')
      connectivityRaw <- many' columnDecimal
      endOfLine
      return
        ( (index, IS.fromList connectivityRaw)
        , Atom
            { _atom_Index        = index
            , _atom_Element      = read cElement
            , _atom_Label        = ""
            , _atom_IsPseudo     = False
            , _atom_FFType       =
                case mFFType of
                  Nothing -> ""
                  Just x' -> show x'
            , _atom_PCharge      = Nothing
            , _atom_Coordinates  = VS.fromList [x, y, z]
            }
          )
    -- Parse multiple non-line-breaking whitespace separated decimals.
    columnDecimal :: Parser Int
    columnDecimal = do
      _ <- many' (char ' ' <|> char '\t')
      i <- decimal
      _ <- many' (char ' ' <|> char '\t')
      return i
-}

{-|
Parse the "interesting" fields of a MOL2 file. This contains partial charges as well as
connectivity. There is no special understanding for the atom types, that are available in MOL2
files. They will simply be treated as the force field string.
-}
parseMOL2 :: Parser Molecule
parseMOL2 = undefined {- do
  (label, nAtoms, nBonds) <- moleculeParser
  atoms <- atomParser
  bonds <- bondParser nAtoms
  let updatedAtoms =
        [ (atoms !! i) & atom_Connectivity .~ I.fromList (bonds !! i)
        | i <- [ 0 .. length atoms - 1 ]
        ]
  return Molecule
    { _molecule_Label    = label
    , _molecule_Atoms    = R.fromList (R.Z R.:. (length updatedAtoms :: Int)) updatedAtoms
    , _molecule_Energy   = Nothing
    , _molecule_Gradient = Nothing
    , _molecule_Hessian  = Nothing
    }
  where
    moleculeParser :: Parser (String, Int, Int)
    moleculeParser = do
      _ <- manyTill anyChar (string "@<TRIPOS>MOLECULE")
      endOfLine
      label <- manyTill anyChar endOfLine
      skipSpace
      nAtoms <- decimal
      _ <- many1 $ char ' '
      nBonds <- decimal
      _ <- manyTill anyChar endOfLine
      return (label, nAtoms, nBonds)
    atomParser :: Parser [Atom]
    atomParser = do
      _ <- manyTill anyChar (string "@<TRIPOS>ATOM")
      endOfLine
      atoms <- many1 atomLineParser
      return atoms
      where
        atomLineParser :: Parser Atom
        atomLineParser = do
          skipSpace
          _ <- decimal
          skipSpace
          label <- manyTill anyChar (char ' ' <|> char '\t')
          skipSpace
          x <- double
          skipSpace
          y <- double
          skipSpace
          z <- double
          skipSpace
          cElement <- many1 letter
          ffDot <- maybeOption $ char '.'
          ffType2 <- maybeOption $ manyTill anyChar (char ' ' <|> char '\t')
          skipSpace
          nSubstructure <- decimal
          skipSpace
          nameSubstructure <- manyTill anyChar (char ' ' <|> char '\t')
          skipSpace
          partialCharge <- double
          _ <- many' (char ' ' <|> char '\t')
          endOfLine
          let ffType =
                if isNothing ffDot || isNothing ffType2
                  then cElement
                  else cElement ++ "." ++ fromJust ffType2
          return Atom
            { _atom_Element      = read cElement
            , _atom_Label        = label
            , _atom_IsPseudo     = False
            , _atom_FFType       = ffType
            , _atom_PCharge      = Just partialCharge
            , _atom_Coordinates  = R.fromListUnboxed (R.Z R.:. 3) [x, y, z]
            , _atom_Connectivity = I.empty
            }
    bondParser :: Int -> Parser [[Int]]
    bondParser nAtoms = do
      _ <- manyTill anyChar (string "@<TRIPOS>BOND")
      endOfLine
      rawBonds <- many1 bondLineParser
      let rawBondsSwaped = map swap rawBonds
          rawBondsComplete = rawBonds ++ rawBondsSwaped
          sortedBonds =
            [ foldr (\(o, t) acc -> if o == i then t : acc else acc) [] rawBondsComplete
            | i <- [ 0 .. nAtoms - 1]
            ]
      return sortedBonds
      where
        bondLineParser :: Parser (Int, Int)
        bondLineParser = do
          skipSpace
          _ <- decimal
          skipSpace
          originAtom <- decimal
          skipSpace
          targetAtom <- decimal
          skipSpace
          _ <- (show <$> decimal) <|> many1 letter
          _ <- many' (char ' ' <|> char '\t')
          endOfLine
          return (originAtom - 1, targetAtom - 1)
-}

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
