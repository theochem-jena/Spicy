module Spicy.Parser
( parse_XYZ
, parse_TXYZ
, parse_MOL2
) where
import           Control.Applicative
import           Data.Attoparsec.Text.Lazy
import           Data.Maybe
import qualified Data.Text                 as T
import           Data.Tuple
import           Lens.Micro.Platform
import           Spicy.Types

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
            Just x  -> show x
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

-- | Parse the "interesting" fields of a MOL2 file. This contains partial
-- | charges as well as connectivity.
parse_MOL2 :: Parser Molecule
parse_MOL2 = do
  (label, nAtoms, nBonds) <- moleculeParser
  atoms <- atomParser
  bonds <- bondParser nAtoms
  let updatedAtoms =
        [ (atoms !! i) & atom_Connectivity .~ (bonds !! i)
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
            , _atom_Connectivity = []
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
