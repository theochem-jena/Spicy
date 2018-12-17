{-
The module collects data only used for unit testing, such as the correct answers
for unit tests.
-}
module Spicy.UnitTests.Data
( moleculeHFeCNxH2O
, moleculeHFeCNxH2OTXYZ
, moleculeHFeCNxH2OXYZ
, moleculeHFeCNxH2OMOL2
, textHFeCNxH2OTXYZ
, testHFeCNxH2OXYZ
, testHFeCNxH2OMOL2
) where
import           Data.Map    (Map)
import qualified Data.Map    as Map
import           Spicy.Types
import qualified Data.Text as T
import Data.Text (Text)
import Lens.Micro.Platform
import Numeric.LinearAlgebra

----------------------------------------------------------------------------------------------------
-- Atom and molecules for testing
----------------------------------------------------------------------------------------------------
-- | Define a molecule HFe(CN)xH2O, with hydrid iron distance being exactly at 1.4 times the sum of
-- | the covalent radii of H and Fe
atomC0 = Atom
  { _atom_Element = C
  , _atom_Label = "C0"
  , _atom_IsPseudo = False
  , _atom_FFType = "6"
  , _atom_PCharge = Just 0.5
  , _atom_Coordinates = ( 0.000960, -0.001240, 0.281420)
  , _atom_Connectivity = [1, 2]
  }

atomN1 = Atom
  { _atom_Element = N
  , _atom_Label = "N1"
  , _atom_IsPseudo = False
  , _atom_FFType = "7"
  , _atom_PCharge = Just (-0.5)
  , _atom_Coordinates = ( 0.000960, -0.001240,-0.875780)
  , _atom_Connectivity = [0]
  }

atomFe2 = Atom
  { _atom_Element = Fe
  , _atom_Label = "Fe2"
  , _atom_IsPseudo = False
  , _atom_FFType = "26"
  , _atom_PCharge = Just 1.4
  , _atom_Coordinates = ( 0.000450, -0.000610, 2.307300)
  , _atom_Connectivity = [0, 3]
  }

atomH3 = Atom
  { _atom_Element = H
  , _atom_Label = "H3"
  , _atom_IsPseudo = False
  , _atom_FFType = "1"
  , _atom_PCharge = Just 0.4
  , _atom_Coordinates = (-0.002370,  0.003010, 4.869300)
  , _atom_Connectivity = [2]
  }

atomO4 = Atom
  { _atom_Element = O
  , _atom_Label = "O4"
  , _atom_IsPseudo = False
  , _atom_FFType = "8"
  , _atom_PCharge = Just (-0.8)
  , _atom_Coordinates = ( 2.059230, -2.830320, 3.280180)
  , _atom_Connectivity = [5, 6]
  }

atomH5 = Atom
  { _atom_Element = H
  , _atom_Label = "H5"
  , _atom_IsPseudo = False
  , _atom_FFType = "1"
  , _atom_PCharge = Just 0.5
  , _atom_Coordinates = ( 2.096560, -2.818240, 2.290710)
  , _atom_Connectivity = [4]
  }

atomH6 = Atom
  { _atom_Element = H
  , _atom_Label = "H6"
  , _atom_IsPseudo = False
  , _atom_FFType = "1"
  , _atom_PCharge = Just 0.4
  , _atom_Coordinates = ( 2.708520, -2.137580, 3.561450)
  , _atom_Connectivity = [4]
  }

-- | Test molecule with all information a molecule could have
moleculeHFeCNxH2O = Molecule
  { _molecule_Label = "HFe(CN)xH2O"
  , _molecule_Atoms =
      [ atomC0
      , atomN1
      , atomFe2
      , atomH3
      , atomO4
      , atomH5
      , atomH6
      ]
  , _molecule_Energy = Just (-1000.0)
  , _molecule_Gradient = Just $ fromList . map fromInteger $ [1 .. 21]
  , _molecule_Hessian = Just $ diagRect 21.0 (fromList [1.0 .. 6.0]) 6 6
  }
-- | Test molecule with all informations a TXYZ could have (without calling Tinker)
moleculeHFeCNxH2OTXYZ = moleculeHFeCNxH2O
  & molecule_Atoms .~ newAtoms
  & molecule_Energy .~ Nothing
  & molecule_Gradient .~ Nothing
  & molecule_Hessian .~ Nothing
  where
    oldAtoms = moleculeHFeCNxH2O ^. molecule_Atoms
    newAtoms = map
      ( (& atom_Label .~ "")
      . (& atom_PCharge .~ Nothing)
      ) oldAtoms

-- | Test molecule with all informations a XYZ could have
moleculeHFeCNxH2OXYZ = moleculeHFeCNxH2O
  & molecule_Atoms .~ newAtoms
  & molecule_Energy .~ Nothing
  & molecule_Gradient .~ Nothing
  & molecule_Hessian .~ Nothing
  where
    oldAtoms = moleculeHFeCNxH2O ^. molecule_Atoms
    newAtoms = map
      ( (& atom_Label .~ "")
      . (& atom_FFType .~ "")
      . (& atom_Connectivity .~ [])
      . (& atom_PCharge .~ Nothing)
      ) oldAtoms

-- | Test molecule with all informations a TXYZ could have
moleculeHFeCNxH2OMOL2 = moleculeHFeCNxH2O
  & molecule_Atoms .~ newAtoms
  & molecule_Energy .~ Nothing
  & molecule_Gradient .~ Nothing
  & molecule_Hessian .~ Nothing
  where
    oldAtoms = moleculeHFeCNxH2O ^. molecule_Atoms
    newFFTypes =
      [ "C.1", "N.pl3", "Fe", "H", "O.3", "H", "H" ]
    newAtoms =
      zipWith (\a ft -> a & atom_FFType .~ ft) oldAtoms newFFTypes

-- | Test molecule as TXYZ file
textHFeCNxH2OTXYZ = T.pack . concat $
  [ "     7 HFe(CN)xH2O\n"
  , "    1  C      0.000960   -0.001240    0.281420     6     2     3     \n"
  , "    2  N      0.000960   -0.001240   -0.875780     7     1           \n"
  , "    3 Fe      0.000450   -0.000610    2.307300    26     1     4     \n"
  , "    4  H     -0.002370    0.003010    4.869300     1     3           \n"
  , "    5  O      2.059230   -2.830320    3.280180     8     6     7     \n"
  , "    6  H      2.096560   -2.818240    2.290710     1     5           \n"
  , "    7  H      2.708520   -2.137580    3.561450     1     5           \n"
  ]

-- | Test molecule as XYZ file
testHFeCNxH2OXYZ = T.pack . concat $
  [ " 7                                         \n"
  , "HFe(CN)xH2O\n"
  , "C         0.000960   -0.001240    0.281420 \n"
  , "N         0.000960   -0.001240   -0.875780 \n"
  , "Fe        0.000450   -0.000610    2.307300 \n"
  , "H        -0.002370    0.003010    4.869300 \n"
  , "O         2.059230   -2.830320    3.280180 \n"
  , "H         2.096560   -2.818240    2.290710 \n"
  , "H         2.708520   -2.137580    3.561450 \n"
  ]

testHFeCNxH2OMOL2 = T.pack . concat $
  [ "@<TRIPOS>MOLECULE\n"
  , "HFe(CN)xH2O\n"
  , " 7 5 0 0 0                                                                    \n"
  , "SMALL                                                                         \n"
  , "GASTEIGER                                                                     \n"
  , "                                                                              \n"
  , "@<TRIPOS>ATOM\n"
  , "      1 C0          0.000960   -0.001240    0.281420 C.1     1  UNL1        0.5000 \n"
  , "      2 N1          0.000960   -0.001240   -0.875780 N.pl3   1  UNL1       -0.5000 \n"
  , "      3 Fe2         0.000450   -0.000610    2.307300 Fe      1  UNL1        1.4000 \n"
  , "      4 H3         -0.002370    0.003010    4.869300 H       1  UNL1        0.4000 \n"
  , "      5 O4          2.059230   -2.830320    3.280180 O.3     1  UNL1       -0.8000 \n"
  , "      6 H5          2.096560   -2.818240    2.290710 H       0  HOH0        0.5000 \n"
  , "      7 H6          2.708520   -2.137580    3.561450 H       0  HOH0        0.4000 \n"
  , "@<TRIPOS>BOND\n"
  , "     1     1     2    1\n"
  , "     2     1     3    3\n"
  , "     3     3     4    1\n"
  , "     4     5     6    1\n"
  , "     5     5     7    1\n"
  , "\n"
  ]
