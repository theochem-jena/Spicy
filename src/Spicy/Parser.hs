module Spicy.Parser
( parseXYZ
, parseTXYZ
, parseMOL2
, parseSpicy
, parseHMatrix
) where
import           Control.Applicative
import           Data.Attoparsec.Text.Lazy
import qualified Data.IntSet               as I
import           Data.Maybe
import qualified Data.Set                  as S
import qualified Data.Text                 as T
import           Data.Tuple
import           Lens.Micro.Platform
import           Spicy.Types
import Numeric.LinearAlgebra hiding (double)
import Data.Char
import Spicy.MolWriter
import qualified Data.Text.IO as T


--------------------------------------------------------------------------------
-- General helper functions
--------------------------------------------------------------------------------
-- | make a parser optional and wrap it in a Maybe
maybeOption :: Parser a -> Parser (Maybe a)
maybeOption p = option Nothing (Just <$> p)


--------------------------------------------------------------------------------
-- Parser for molecular formats
--------------------------------------------------------------------------------
-- | parse a .xyz file (has no connectivity, atom types or partioal charges)
parseXYZ :: Parser Molecule
parseXYZ = do
  skipSpace
  nAtoms <- decimal
  skipSpace
  comment <- manyTill anyChar endOfLine
  atoms <- count nAtoms xyzLineParser
  return Molecule
    { _molecule_Label = comment
    , _molecule_Atoms = atoms
    , _molecule_Energy = Nothing
    , _molecule_Gradient = Nothing
    , _molecule_Hessian = Nothing
    }
  where
    xyzLineParser :: Parser Atom
    xyzLineParser = do
      skipSpace
      cElement <- many1 letter
      skipSpace
      x <- double
      skipSpace
      y <- double
      skipSpace
      z <- double
      skipSpace
      return Atom
        { _atom_Element = read cElement
        , _atom_Label = ""
        , _atom_IsPseudo = False
        , _atom_FFType = ""
        , _atom_PCharge = Nothing
        , _atom_Coordinates = (x, y, z)
        , _atom_Connectivity = I.empty
        }

-- | parse a .txyz file (Tinkers xyz format)
-- | it has coordinates and might have connectivity and atom types
parseTXYZ :: Parser Molecule
parseTXYZ = do
  skipSpace
  nAtoms <- decimal
  _ <- many' (char ' ' <|> char '\t')
  comment <- manyTill anyChar endOfLine
  atoms <- many1 txyzLineParser
  return Molecule
    { _molecule_Label = comment
    , _molecule_Atoms = atoms
    , _molecule_Energy = Nothing
    , _molecule_Gradient = Nothing
    , _molecule_Hessian = Nothing
    }
  where
    txyzLineParser :: Parser Atom
    txyzLineParser = do
      skipSpace
      _ <- decimal
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
      return Atom
        { _atom_Element = read cElement
        , _atom_Label = ""
        , _atom_IsPseudo = False
        , _atom_FFType =
          case mFFType of
            Nothing -> ""
            Just x  -> show x
        , _atom_PCharge = Nothing
        , _atom_Coordinates = (x, y, z)
        , _atom_Connectivity = I.fromList $ map (+ (-1)) connectivityRaw
        }
    columnDecimal :: Parser Int
    columnDecimal = do
      _ <- many' (char ' ' <|> char '\t')
      i <- decimal
      _ <- many' (char ' ' <|> char '\t')
      return i

-- | Parse the "interesting" fields of a MOL2 file. This contains partial
-- | charges as well as connectivity.
parseMOL2 :: Parser Molecule
parseMOL2 = do
  (label, nAtoms, nBonds) <- moleculeParser
  atoms <- atomParser
  bonds <- bondParser nAtoms
  let updatedAtoms =
        [ (atoms !! i) & atom_Connectivity .~ I.fromList (bonds !! i)
        | i <- [ 0 .. length atoms - 1 ]
        ]
  return Molecule
    { _molecule_Label = label
    , _molecule_Atoms = updatedAtoms
    , _molecule_Energy = Nothing
    , _molecule_Gradient = Nothing
    , _molecule_Hessian = Nothing
    }
  where
    moleculeParser :: Parser (String, Int, Int)
    moleculeParser = do
      _ <- manyTill anyChar (string $ T.pack "@<TRIPOS>MOLECULE")
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
      _ <- manyTill anyChar (string $ T.pack "@<TRIPOS>ATOM")
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
            { _atom_Element = read cElement
            , _atom_Label = label
            , _atom_IsPseudo = False
            , _atom_FFType = ffType
            , _atom_PCharge = Just partialCharge
            , _atom_Coordinates = (x, y, z)
            , _atom_Connectivity = I.empty
            }
    bondParser :: Int -> Parser [[Int]]
    bondParser nAtoms = do
      _ <- manyTill anyChar (string $ T.pack "@<TRIPOS>BOND")
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

-- | Parser for the Spicy format used in this program
parseSpicy :: Parser Molecule
parseSpicy = do
  _ <- string $ T.pack "#Spicy-Format v0.1"
  endOfLine
  skipSpace
  _ <- string $ T.pack "#Spicy-Molecule"
  endOfLine
  _ <- string $ T.pack "  Label:"
  endOfLine
  _ <- string $ T.pack "    "
  mLabel <- takeTill isEndOfLine
  skipSpace
  mEnergy <- maybeOption parseEnergy
  mGradient <- maybeOption parseGradient
  mHessian <- maybeOption parseHessian
  skipSpace
  _ <- string $ T.pack "#Spicy-Atoms"
  skipSpace
  mAtoms <- many1 parseAtoms
  return Molecule
    { _molecule_Label = T.unpack mLabel
    , _molecule_Atoms = mAtoms
    , _molecule_Energy = mEnergy
    , _molecule_Gradient = mGradient --  mGradient
    , _molecule_Hessian = mHessian -- mHessian
    }
  where
    parseEnergy = do
      _ <- manyTill anyChar (string $ T.pack "Energy / Hartree:")
      skipSpace
      energy <- double
      return energy
    parseGradient = do
      _ <- manyTill anyChar (string $ T.pack "Gradient / Hartee/Bohr:")
      skipSpace
      gradient <- many1 $ do
        skipSpace
        gVal <- double
        return gVal
      return $ fromList gradient
    parseHessian = do
      _ <- manyTill anyChar (string $ T.pack "Hessian / a.u.:")
      skipSpace
      hessian <- parseHMatrix
      return hessian
    parseAtoms = do
      _ <- many' (char ' ' <|> char '\t')
      cElement <- many1 letter
      _ <- many1 (char ' ' <|> char '\t')
      label <- takeTill isSpace
      _ <- many1 (char ' ' <|> char '\t')
      pseudo <- option ' ' (char 'P')
      _ <- many' (char ' ' <|> char '\t')
      ffType <- takeTill isHorizontalSpace
      _ <- many1 (char ' ' <|> char '\t')
      pChargeTest <- maybeOption $ string $ T.pack "No"
      pCharge <- if pChargeTest == Nothing
        then Just <$> double
        else return Nothing
      _ <- many1 (char ' ' <|> char '\t')
      coordVec <- count 3 $ do
        coordComponent <- double
        _ <- many1 (char ' ' <|> char '\t')
        return coordComponent
      _ <- many' (char ' ' <|> char '\t')
      connectivity <- many' $ do
        conAtom <- decimal
        _ <- many' (char ' ' <|> char '\t')
        return conAtom
      endOfLine
      return Atom
        { _atom_Element = read cElement
        , _atom_Label = T.unpack label
        , _atom_IsPseudo = if pseudo == 'P' then True else False
        , _atom_FFType = T.unpack ffType
        , _atom_PCharge = pCharge
        , _atom_Coordinates = (\[x, y, z] -> (x, y, z)) coordVec
        , _atom_Connectivity = I.fromList connectivity
        }

main = do
  raw <- T.readFile "/home/phillip/Spicy.spc"
  let spc = parseOnly parseSpicy raw
  putStrLn $ show spc

----------------------------------------------------------------------------------------------------
-- Parser for generic formats
----------------------------------------------------------------------------------------------------
-- | Parse the "show" instance output for HMatrix' matrix Type, including the dimension infos
parseHMatrix :: Parser (Matrix Double)
parseHMatrix = do
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
  return $ fromRows . map fromList $ rows
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
