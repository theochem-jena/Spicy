{-|
Module      : Spicy.MolWriter
Description : Writing chemical data formats
Copyright   : Phillip Seeber, 2019
License     : GPL-3
Maintainer  : phillip.seeber@uni-jena.de
Stability   : experimental
Portability : POSIX, Windows

A module which converts the internal Molecule representation to a string, which is a common chemical
file format, that can be read by Avogadro, VMD, OpenBabel etc.. The writers are not fool proof with
respect to force field types, which should always be remembered when usings its results.
-}
{-# LANGUAGE ScopedTypeVariables #-}
module Spicy.MolWriter
( writeXYZ
, writeTXYZ
, writeMOL2
, writePDB
, writeSpicy
) where
import           Data.Aeson.Encode.Pretty
import           Data.Attoparsec.Text.Lazy (isEndOfLine)
import qualified Data.IntMap.Lazy          as IM
import           Data.Maybe
import qualified Data.Sequence             as S
import           Data.Text.Lazy            (Text)
import qualified Data.Text.Lazy            as T
import qualified Data.Text.Lazy.Encoding   as T
import           Lens.Micro.Platform
import           Prelude                   hiding (cycle, foldl1, foldr1, head,
                                            init, last, maximum, minimum, tail,
                                            take, takeWhile, (!!))
import           Spicy.Types
import           Text.Printf
import Data.Foldable


{-|
Write a Molden XYZ file from a molecule. This format ignores all deep level layers of a molecule.
-}
writeXYZ :: Molecule -> Either String Text
writeXYZ mol
  | coordCheck = Right xyzText
  | otherwise  =
      Left "writeXYZ: Atomic coordinates damaged. Some of the atoms seem to have not 3 coordinates."
  where
    -- Check if each atom has exactly 3 coordinates.
    coordCheck = all (\a -> S.length (a ^. atom_Coordinates) == 3) $ mol ^. molecule_Atoms
    -- Line 1 of the header contains the number of atoms
    headerL1   = T.pack . show . IM.size $ mol ^. molecule_Atoms
    -- Line 2 of the header contains 1 line of comment. Remove all linebreaks by filtering.
    headerL2   = T.filter (not . isEndOfLine) $ mol ^. molecule_Label
    atomLs     =
      IM.foldr' (\atom acc ->
        (T.pack $ printf "%-4s    %12.8F    %12.8F    %12.8F\n"
          (show $ atom ^. atom_Element)
          (fromMaybe 0.0 $ (atom ^. atom_Coordinates) S.!? 0)
          (fromMaybe 0.0 $ (atom ^. atom_Coordinates) S.!? 1)
          (fromMaybe 0.0 $ (atom ^. atom_Coordinates) S.!? 2)
        )
        `T.append`
        acc
      ) "" $ mol ^. molecule_Atoms
    -- The assembled text, repesenting a valid XYZ file.
    xyzText    =
      T.unlines
        [ headerL1
        , headerL2
        ]
      `T.append`
      atomLs

{-|
Write a Tinker XYZ from a 'Molecule'. The writer trusts the '_atom_FFType' to be
correct (if set) and will simply write them out. Therefore it is possible, that wrong atom types can
be written. If they are not set, the writer will simply equalise all atom types to 0, which is OK
for visualisation but obviously not for MM.

This format ingores all deeper level layers of a molecule.
-}
writeTXYZ :: Molecule -> Text
writeTXYZ _mol = undefined
{-
  let atoms        = mol ^. molecule_Atoms
      nAtoms       = VB.length atoms
      -- Index the atoms. Ignore higher
      atomsIndexed = VB.zip (VB.generate nAtoms (\i -> i)) atoms
  in  -- Header line with number of atoms and comment
      ( (T.pack $ show nAtoms ++ "  " ++ (mol ^. molecule_Label))
        `T.append`
         nL
      )
      `T.append`
      -- Body with one atom per line
      ( vUnlines . VB.map (\(i, a) ->
          (T.pack . printf "%-6d  " $ i + 1)
          `T.append`
          a2txyz a
          `T.append`
          m2txyz mol i
        ) $ atomsIndexed
      )
  where
    -- Write element, XYZ coordinates and force field type of an atom.
    a2txyz :: Atom -> Text
    a2txyz a =
      (T.pack . printf "%-4s" . show $ a ^. atom_Element)
      `T.append`
      (vConcat . VB.map (T.pack . printf "    %12.8F") . VS.convert $ a ^. atom_Coordinates)
      `T.append`
      "      "
      `T.append`
      ( if a ^. atom_FFType == ""
          then T.pack . printf "%10s     " $ ("0" :: String)
          else T.pack . printf "%10s     " $ a ^. atom_FFType
      )
    -- Write connectivity of an atom (specified by its index). If no bonding informations can be
    -- found, do not write them.
    m2txyz :: Molecule -> Int -> Text
    m2txyz m i =
      let iBonds = i `IM.lookup` (m ^. molecule_Bonds)
      in  case iBonds of
            Just b  ->
              IS.foldl' (\acc x -> acc `T.append` (T.pack $ printf "%6d  " x)) "" (IS.map (+1) b)
            Nothing -> ""
-}

{-|
Write a simplified .mol2 file (Tripos SYBYL) from a 'Molecule', containing the atoms, connectivities
(single bonds only) and partial charges. The writer is not fool proof and will happily accept any
'_atom_FFType' that is supplied, even if it is not a TRIPOS SYBYL atom type. This can lead to mol2
files that have correct topology and geometry, but visualisation programs wont be able to assign
correct elements.
-}
writeMOL2 :: Molecule -> Text
writeMOL2 mol = T.decodeUtf8 . encodePretty $ mol
{-
  let atoms  = mol ^. molecule_Atoms
      nAtoms = VB.length atoms
      bonds  = mol ^. molecule_Bonds
      nBonds = IM.size bonds
  in  -- Header of the MOL2.
      T.unlines
        [ "@<TRIPOS>MOLECULE"
        , T.pack $ mol ^. molecule_Label
        , T.pack $ show nAtoms ++ " " ++ show nBonds ++ " 0 0 0"
        , "SMALL"
        , "GASTEIGER"
        , ""
        ]
      `T.append`
      -- Atoms section containing indices, chemical element, XYZ coordinates, force field type
      "@<TRIPOS>ATOM"
      `T.append`
  where
    a2mol :: Atom -> Text
    a2mol a
-}


{-
  "@<TRIPOS>MOLECULE" ++ "\n" ++
  mol ^. molecule_Label ++ "\n" ++
  show nAtoms ++ " " ++ show nBonds ++ " 0 0 0" ++ "\n" ++
  "SMALL" ++ "\n" ++
  "GASTEIGER" ++ "\n\n" ++

  "@<TRIPOS>ATOM" ++ "\n" ++
  concat
    ( map (\(n, a) ->
        printf "%6d    "        n ++
        printf "%-4s    "       (show $ a ^. atom_Element) ++
        (\(x, y, z) -> printf "%12.8F    %12.8F    %12.8F        " x y z)
          (indexAtomCoordinates $ a ^. atom_Coordinates) ++
        printf "%-8s    "
          (if a ^. atom_FFType == ""
             then show (a ^. atom_Element) ++
                  "." ++
                  show (length . I.toList $ a ^. atom_Connectivity)
             else a ^. atom_FFType
          ) ++
        printf "%2d    " (1 :: Int) ++
        printf "%4s    " "UNL1" ++
        printf "%12.8F\n" (fromMaybe 0.0 $ a ^. atom_PCharge)
      ) numberedAtoms
    ) ++ "\n" ++
    "@<TRIPOS>BOND" ++ "\n" ++
    concat
      ( map (\(n, (o, t)) -> printf "    %6d    " n ++
                           printf "%6d    " o ++
                           printf "%6d    " t ++
                           printf "%6d\n" (1 :: Int)
            ) numberedBonds
      )
  where
    atoms = mol ^. molecule_Atoms
    nAtoms = V.length atoms
    numberedAtoms = V.generate nAtoms (\i -> (i, atoms V.! i))
    bonds = map (I.toList . (^. atom_Connectivity)) $ V.toList atoms
    pairBondsRedundant =
      concat
      [ map (\a -> (i + 1, a + 1)) (bonds !! i)
      | i <- [ 0 .. length bonds - 1 ]
      ]
    pairBonds =
      foldr (\a acc -> if swap a `elem` acc
                       then delete a acc
                       else acc
            ) pairBondsRedundant pairBondsRedundant
    nBonds = length pairBonds
    bondIndexList = [ 1 .. nBonds]
    numberedBonds = V.generate nBonds (\i -> (i, pairBonds !! i))
-}

writePDB :: Molecule -> Text
writePDB _mol = undefined

{-|
Write Spicy format, which is an AESON generated JSON document, directly representing the data type
of 'Molecule'.
-}
writeSpicy :: Molecule -> Text
writeSpicy mol = T.decodeUtf8 . encodePretty $ mol
