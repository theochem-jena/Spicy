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
{-# LANGUAGE OverloadedStrings #-}
module Spicy.MolWriter
( writeXYZ
, writeTXYZ
, writeMOL2
, writeSpicy
, Molecule(..)
) where
import qualified Data.Array.Accelerate             as A
import qualified Data.Array.Accelerate.LLVM.Native as A
import qualified Data.Array.Accelerate.IO.Data.Vector.Generic as A
import qualified Data.IntSet                       as I
import           Data.List
import           Data.List.Split
import           Data.Maybe
import           Data.Text.Lazy                    (Text)
import qualified Data.Text.Lazy                    as T
import           Data.Tuple
import qualified Data.Vector                       as VB
import qualified Data.Vector.Storable as VS
import           Lens.Micro.Platform
import           Spicy.Types
import           Text.Printf

{-|
Portable new line character.
-}
nL :: Text
nL = T.unlines [""]

{-|
Write a .xyz file from a molecule.
-}
writeXYZ :: Molecule -> Text
writeXYZ m =
  -- Header section with number of atoms and comment
  T.unlines (
    [ T.pack . show $ (VB.length $ m ^. molecule_Atoms)
    , T.pack $ m ^. molecule_Label
    ]
  )
  `T.append`
  -- Body with element symbols and XYZ coordinates
  (VB.foldl' (\acc x -> acc `T.append` x `T.append` nL) "" $ VB.map a2xyz $ m ^. molecule_Atoms)
  where
    -- Write informations about a single atom to a line
    a2xyz :: Atom -> Text
    a2xyz a =
      (T.pack . printf "%-4s" . show $ a ^. atom_Element)
      `T.append`
      VB.foldl' (T.append) "" (VB.map (T.pack . printf "    %12.8F") . VS.convert $ a ^. atom_Coordinates)
{-
  show nAtoms ++ "\n" ++
  comment ++ "\n" ++
  concat
    ( map ( \a ->
        printf "%-4s    "   (show (a ^. atom_Element)) ++
        (\(x, y, z) -> printf "%12.8F    %12.8F    %12.8F\n" x y z)
          (indexAtomCoordinates $ a ^. atom_Coordinates)
        {-
        printf "%12.8F    " ((a ^. atom_Coordinates) ! (A.Z A.:. 0)) ++
        printf "%12.8F    " ((a ^. atom_Coordinates) ! (A.Z A.:. 1)) ++
        printf "%12.8F\n"   ((a ^. atom_Coordinates) ! (A.Z A.:. 2))
        -}
      ) $ V.toList (mol ^.  molecule_Atoms)
    )
  where
    nAtoms = V.length $ mol ^. molecule_Atoms
    comment = mol ^. molecule_Label
-}

{-|
Write a .txyz (Tinkers xyz format) from a 'Molecule'. The writer trusts the '_atom_FFType' to be
correct (if set) and will simply write them out. Therefore it is possible, that wrong atom types can
be written. If they are not set, the writer will simply equalise all atom types to 0, which is OK
for visualisation but obviously not for MM.
-}
writeTXYZ :: Molecule -> String
writeTXYZ mol = undefined
{-
  show nAtoms ++ "  " ++ comment ++ "\n" ++
  concat
    ( map (\(n, a) ->
        printf "%-6d  "         n ++
        printf "%-4s    "       (show (a ^. atom_Element)) ++
        (\(x, y, z) -> printf "%12.8F    %12.8F    %12.8F        " x y z)
          (indexAtomCoordinates $ a ^. atom_Coordinates) ++
        printf "%6s      "
          (if a ^. atom_FFType == ""
            then "0"
            else a ^. atom_FFType
          ) ++
        concat
          ( map (printf "%6d  " . (+ 1)) (I.toList $ a ^. atom_Connectivity)
          ) ++ "\n"
      ) $ V.toList numberedAtoms
    )
  where
    nAtoms = V.length $ mol ^. molecule_Atoms
    comment = head . lines $ mol ^. molecule_Label
    atoms = mol ^. molecule_Atoms
    numberedAtoms = V.generate nAtoms (\i -> (i, atoms V.! i))
-}

{-|
Write a simplified .mol2 file (Tripos SYBYL) from a 'Molecule', containing the atoms, connectivities
(single bonds only) and partial charges. The writer is not fool proof and will happily accept any
'_atom_FFType' that is supplied, even if it is not a TRIPOS SYBYL atom type. This can lead to mol2
files that have correct topology and geometry, but visualisation programs wont be able to assign
correct elements.
-}
writeMOL2 :: Molecule -> String
writeMOL2 mol = undefined
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

{-|
Write Spicy format, which is a custom format (not stable yet), containing all informations, that are
internally used to represent a molecule.
-}
writeSpicy :: Molecule -> String
writeSpicy m = undefined
{-
  "#Spicy-Format v0.2\n" ++
  "\n" ++
  "#Spicy-Molecule\n" ++
  "  Label:\n" ++
  "    " ++ (m ^. molecule_Label) ++ "\n" ++
  case energy of
    Nothing -> ""
    Just e ->
      "  Energy / Hartree:\n" ++
      "    " ++ printf "%20.10e\n" e
  ++
  case gradient of
    Nothing -> ""
    Just g ->
      "  Gradient / Hartee/Bohr:\n" ++
      (concat . map ("  " ++) . map writeSafeListToLine . chunksOf 3 . V.toList $ g)
      --( concat . map ((++ "\n") . ("    " ++) . show) . chunksOf 3 . toList $ g )
  ++
  case hessian of
    Nothing -> ""
    Just h ->
      "  Hessian / a.u.:\n" ++
      ( concat . map ((++ "\n") . ("    " ++)) . splitOn "\n" . show $ h )
  ++
  "\n" ++
  "#Spicy-Atoms\n" ++
  ( concat . map
      ( \a ->
          printf "  %3s  " (show $ a ^. atom_Element) ++
          printf "  %6s  " (a ^. atom_Label) ++
          printf "  %1s  " (if (a ^. atom_IsPseudo) then "P" else "") ++
          printf "  %6s  " (a ^. atom_FFType) ++
          printf "  %8s  "
            ( case (a ^. atom_PCharge) of
                Nothing -> "No"
                Just c  -> printf "%8.5f" c
            ) ++
          (\(x, y, z) -> printf "  %16.10f  %16.10f  %16.10f  " x y z)
            (indexAtomCoordinates $ a ^. atom_Coordinates) ++
          (concat . map (printf "  %6d  ") . I.toList $ a ^. atom_Connectivity) ++
          "\n"
      ) $ V.toList atoms
    )
  where
    atoms = m ^. molecule_Atoms
    energy = m ^. molecule_Energy
    gradient = m ^. molecule_Gradient
    hessian = m ^. molecule_Hessian
    writeSafeListToLine :: (Floating a, PrintfArg a) => [a] -> String
    writeSafeListToLine [] = ""
    writeSafeListToLine [x] = printf "  %8.5f\n" x
    writeSafeListToLine (x:xs) =
      printf "  %8.5f  " x ++
      writeSafeListToLine xs
-}
