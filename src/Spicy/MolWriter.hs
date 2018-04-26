module Spicy.MolWriter
( write_XYZ
, write_TXYZ
) where
import           Lens.Micro.Platform
import           Spicy.Types
import           Text.Printf

-- | writing a .xyz file from a molecule
write_XYZ :: Molecule -> String
write_XYZ mol =
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

-- | write a .txyz (Tinkers xyz format) from a molecule
write_TXYZ :: Molecule -> String
write_TXYZ mol =
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
