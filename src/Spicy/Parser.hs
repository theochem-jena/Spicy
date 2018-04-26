module Spicy.Parser
( parse_XYZ
, parse_TXYZ
) where
import Data.Attoparsec.Text.Lazy
import Spicy.Types
import Control.Applicative


-- | make a parser optional
maybeOption :: Parser a -> Parser (Maybe a)
maybeOption p = option Nothing (Just <$> p)


--------------------------------------------------------------------------------
-- parser for molecular formats
--------------------------------------------------------------------------------
-- | parse a .xyz file (has no connectivity, atom types or partioal charges)
parse_XYZ :: Parser Molecule
parse_XYZ = do
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
        , _atom_Connectivity = []
        }

-- | parse a .txyz file (Tinkers xyz format)
-- | it has coordinates and might have connectivity and atom types
parse_TXYZ :: Parser Molecule
parse_TXYZ = do
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
            Just x -> show x
        , _atom_PCharge = Nothing
        , _atom_Coordinates = (x, y, z)
        , _atom_Connectivity = map (+ (-1)) connectivityRaw
        }
    columnDecimal :: Parser Int
    columnDecimal = do
      _ <- many' (char ' ' <|> char '\t')
      i <- decimal
      _ <- many' (char ' ' <|> char '\t')
      return i
