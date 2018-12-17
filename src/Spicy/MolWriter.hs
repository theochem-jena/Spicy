{-
A module which converts the internal Molecule representation to a string, which
is a common chemical file format, that can be read by Avogadro, VMD, OpenBabel
etc.. The writers are not fool proof with respect to force field types, which
should always be remembered when usings its results.
-}
module Spicy.MolWriter
( writeXYZ
, writeTXYZ
, writeMOL2
, writeSpicy
) where
import           Data.List
import           Data.Maybe
import           Data.Tuple
import           Lens.Micro.Platform
import           Spicy.Types
import           Text.Printf
import Numeric.LinearAlgebra
import Data.List.Split


--------------------------------------------------------------------------------
-- Converters from Molecules to chemical formats
--------------------------------------------------------------------------------
-- | Write a .xyz file from a molecule
writeXYZ :: Molecule -> String
writeXYZ mol =
  show nAtoms ++ "\n" ++
  comment ++ "\n" ++
  concat
    ( map ( \a -> printf "%-4s    "   (show (a ^. atom_Element)) ++
                  printf "%12.8F    " (a ^. atom_Coordinates . _1) ++
                  printf "%12.8F    " (a ^. atom_Coordinates . _2) ++
                  printf "%12.8F\n"   (a ^. atom_Coordinates . _3)
          ) (mol ^. molecule_Atoms)
    )
  where
    nAtoms = length (mol ^. molecule_Atoms)
    comment = head . lines $ mol ^. molecule_Label

-- | Write a .txyz (Tinkers xyz format) from a molecule. The writer trusts the
-- | _atom_FFType to be correct (if set) and will simply write them out.
-- | Therefore it is possible, that wrong atom types can be written. If they are
-- | not set, the writer will simply equalize all atom types to 0, which is OK
-- | for visualisation but obviously not for MM
writeTXYZ :: Molecule -> String
writeTXYZ mol =
  show nAtoms ++ "  " ++ comment ++ "\n" ++
  concat
    ( map (\(n, a) -> printf "%-6d  "         n ++
                      printf "%-4s    "       (show (a ^. atom_Element)) ++
                      printf "%12.8F    "     (a ^. atom_Coordinates . _1) ++
                      printf "%12.8F    "     (a ^. atom_Coordinates . _2) ++
                      printf "%12.8F        " (a ^. atom_Coordinates . _3) ++
                      printf "%6s      "
                        (if a ^. atom_FFType == ""
                          then "0"
                          else a ^. atom_FFType
                        ) ++
                      concat
                        ( map (printf "%6d  " . (+ 1)) (a ^. atom_Connectivity)
                        ) ++ "\n"
          ) numberedAtoms
    )

  where
    nAtoms = length (mol ^. molecule_Atoms)
    comment = head . lines $ mol ^. molecule_Label
    atoms = mol ^. molecule_Atoms
    atomIndexList = [ 1 .. length atoms ]
    numberedAtoms = zip atomIndexList atoms

-- | Write a simplified .mol2 file (Tripos SYBYL) from a molecule, containing
-- | the atoms, connectivities (single bonds only) and partial charges.
-- | The writer is not fool proof and will happily accept any _atom_FFType that
-- | is supplied, even if it not even TRIPOS SYBYL atom type. This can lead to
-- | mol2 files that have correct topology and geometry, but visualisation
-- | programs wont be able to assign correct elements.
writeMOL2 :: Molecule -> String
writeMOL2 mol =
  "@<TRIPOS>MOLECULE" ++ "\n" ++
  mol ^. molecule_Label ++ "\n" ++
  show nAtoms ++ " " ++ show nBonds ++ " 0 0 0" ++ "\n" ++
  "SMALL" ++ "\n" ++
  "GASTEIGER" ++ "\n\n" ++

  "@<TRIPOS>ATOM" ++ "\n" ++
  concat
    ( map (\(n, a) -> printf "%6d    " n ++
                      printf "%-4s    " (show $ a ^. atom_Element) ++
                      printf "%12.8F    " (a ^. atom_Coordinates . _1) ++
                      printf "%12.8F    " (a ^. atom_Coordinates . _2) ++
                      printf "%12.8F        " (a ^. atom_Coordinates . _3) ++
                      printf "%-8s    "
                        (if a ^. atom_FFType == ""
                           then show (a ^. atom_Element) ++
                                "." ++
                                show (length $ a ^. atom_Connectivity)
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
    atomIndexList = [ 1 .. nAtoms ]
    numberedAtoms = zip atomIndexList atoms
    nAtoms = length atoms
    bonds = map (^. atom_Connectivity) atoms
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
    numberedBonds = zip bondIndexList pairBonds

-- | Write Spicy format, which is a custom format (not stable yet), containing
-- | all informations, that are internally used to represent a molecule
writeSpicy :: Molecule -> String
writeSpicy m =
  "#Spicy-Format v0.1\n" ++
  "\n" ++
  "#Spicy-Molecule\n" ++
  "  Label:\n" ++
  "    " ++ show (m ^. molecule_Label) ++ "\n" ++
  case energy of
    Nothing -> ""
    Just e ->
      "  Energy / Hartree:\n" ++
      "    " ++ printf "%20.10E\n" (show $ m ^. molecule_Energy)
  ++
  case gradient of
    Nothing -> ""
    Just g ->
      "  Gradient / Hartee/Bohr:\n" ++
      ( concat . map ((++ "\n") . ("    " ++) . show) . chunksOf 3 . toList $ g )
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
                Just c -> printf "%8.5f" c
            ) ++
          printf "  %16.10F  %16.10F  %16.10F  "
            (a ^. atom_Coordinates ^. _1)
            (a ^. atom_Coordinates ^. _2)
            (a ^. atom_Coordinates ^. _3) ++
          (concat . map (printf "  %6d  ") $ a ^. atom_Connectivity) ++
          "\n"
      ) $ atoms
    )
  where
    atoms = m ^. molecule_Atoms
    energy = m ^. molecule_Energy
    gradient = m ^. molecule_Gradient
    hessian = m ^. molecule_Hessian
