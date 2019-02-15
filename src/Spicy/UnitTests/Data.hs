{-
The module collects data only used for unit testing, such as the correct answers
for unit tests.
-}
module Spicy.UnitTests.Data
( moleculeEmpty
, moleculeHFeCNxH2O
, moleculeHFeCNxH2OTXYZ
, moleculeHFeCNxH2OXYZ
, moleculeHFeCNxH2OMOL2
, textHFeCNxH2OTXYZ
, testHFeCNxH2OXYZ
, testHFeCNxH2OMOL2
, testHFeCNxH2OSpicy
) where
import qualified Data.IntSet           as I
import           Data.Text             (Text)
import qualified Data.Text             as T
import           Lens.Micro.Platform
import           Numeric.LinearAlgebra
import           Spicy.Types

----------------------------------------------------------------------------------------------------
-- Atom and molecules for testing
----------------------------------------------------------------------------------------------------
-- | Empty molecule
moleculeEmpty :: Molecule
moleculeEmpty = Molecule
  { _molecule_Label = ""
  , _molecule_Atoms = []
  , _molecule_Energy = Nothing
  , _molecule_Gradient = Nothing
  , _molecule_Hessian = Nothing
  }

-- | Define a molecule HFe(CN)xH2O, with hydrid iron distance being exactly at 1.4 times the sum of
-- | the covalent radii of H and Fe
atomC0 :: Atom
atomC0 = Atom
  { _atom_Element = C
  , _atom_Label = "C0"
  , _atom_IsPseudo = False
  , _atom_FFType = "6"
  , _atom_PCharge = Just 0.5
  , _atom_Coordinates = ( 0.000960, -0.001240, 0.281420)
  , _atom_Connectivity = I.fromList [1, 2]
  }

atomN1 :: Atom
atomN1 = Atom
  { _atom_Element = N
  , _atom_Label = "N1"
  , _atom_IsPseudo = False
  , _atom_FFType = "7"
  , _atom_PCharge = Just (-0.5)
  , _atom_Coordinates = ( 0.000960, -0.001240,-0.875780)
  , _atom_Connectivity = I.fromList [0]
  }

atomFe2 :: Atom
atomFe2 = Atom
  { _atom_Element = Fe
  , _atom_Label = "Fe2"
  , _atom_IsPseudo = False
  , _atom_FFType = "26"
  , _atom_PCharge = Just 1.4
  , _atom_Coordinates = ( 0.000450, -0.000610, 2.307300)
  , _atom_Connectivity = I.fromList [0, 3]
  }

atomH3 :: Atom
atomH3 = Atom
  { _atom_Element = H
  , _atom_Label = "H3"
  , _atom_IsPseudo = True
  , _atom_FFType = "1"
  , _atom_PCharge = Just 0.4
  , _atom_Coordinates = (-0.002370,  0.003010, 4.869300)
  , _atom_Connectivity = I.fromList [2]
  }

atomO4 :: Atom
atomO4 = Atom
  { _atom_Element = O
  , _atom_Label = "O4"
  , _atom_IsPseudo = False
  , _atom_FFType = "8"
  , _atom_PCharge = Just (-0.8)
  , _atom_Coordinates = ( 2.059230, -2.830320, 3.280180)
  , _atom_Connectivity = I.fromList [5]
  }

atomH5 :: Atom
atomH5 = Atom
  { _atom_Element = H
  , _atom_Label = "H5"
  , _atom_IsPseudo = False
  , _atom_FFType = "1"
  , _atom_PCharge = Nothing
  , _atom_Coordinates = ( 2.096560, -2.818240, 2.290710)
  , _atom_Connectivity = I.fromList [4]
  }

atomH6 :: Atom
atomH6 = Atom
  { _atom_Element = H
  , _atom_Label = "H6"
  , _atom_IsPseudo = True
  , _atom_FFType = "1"
  , _atom_PCharge = Just 0.4
  , _atom_Coordinates = ( 2.708520, -2.137580, 3.561450)
  , _atom_Connectivity = I.empty
  }

-- | Test molecule with all information a molecule could have
moleculeHFeCNxH2O :: Molecule
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
moleculeHFeCNxH2OTXYZ :: Molecule
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
      . (& atom_IsPseudo .~ False)
      ) oldAtoms

-- | Test molecule with all informations a XYZ could have
moleculeHFeCNxH2OXYZ :: Molecule
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
      . (& atom_Connectivity .~ I.empty)
      . (& atom_PCharge .~ Nothing)
      . (& atom_IsPseudo .~ False)
      ) oldAtoms

-- | Test molecule with all informations a TXYZ could have
moleculeHFeCNxH2OMOL2 :: Molecule
moleculeHFeCNxH2OMOL2 = moleculeHFeCNxH2O
  & molecule_Atoms .~ newAtoms
  & molecule_Energy .~ Nothing
  & molecule_Gradient .~ Nothing
  & molecule_Hessian .~ Nothing
  where
    oldAtoms = moleculeHFeCNxH2O ^. molecule_Atoms
    newFFTypes =
      [ "C.1", "N.pl3", "Fe", "H", "O.3", "H", "H" ]
    newAtoms'' =
      zipWith (\a ft -> a & atom_FFType .~ ft) oldAtoms newFFTypes
    newAtoms' = map
      ( (& atom_IsPseudo .~ False)
      ) newAtoms''
    newAtoms = newAtoms' & (ix 5 . atom_PCharge) .~ Just 0.5


textHFeCNxH2OTXYZ :: Text-- | Test molecule as TXYZ file
textHFeCNxH2OTXYZ = T.pack . concat $
  [ "     7 HFe(CN)xH2O\n"
  , "    1  C      0.000960   -0.001240    0.281420     6     2     3     \n"
  , "    2  N      0.000960   -0.001240   -0.875780     7     1           \n"
  , "    3 Fe      0.000450   -0.000610    2.307300    26     1     4     \n"
  , "    4  H     -0.002370    0.003010    4.869300     1     3           \n"
  , "    5  O      2.059230   -2.830320    3.280180     8     6           \n"
  , "    6  H      2.096560   -2.818240    2.290710     1     5           \n"
  , "    7  H      2.708520   -2.137580    3.561450     1                 \n"
  ]


testHFeCNxH2OXYZ :: Text-- | Test molecule as XYZ file
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

testHFeCNxH2OMOL2 :: Text
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
  , "         1         1         2         1\n"
  , "         2         1         3         1\n"
  , "         3         3         4         1\n"
  , "         4         5         6         1\n"
  , "\n"
  ]

testHFeCNxH2OSpicy :: Text
testHFeCNxH2OSpicy = T.pack . concat $
  [ "#Spicy-Format v0.2\n"
  , "\n"
  , "#Spicy-Molecule\n"
  , "  Label:\n"
  , "    HFe(CN)xH2O\n"
  , "  Energy / Hartree:\n"
  , "         -1.0000000000e3\n"
  , "  Gradient / Hartee/Bohr:\n"
  , "     1.00000     2.00000     3.00000\n"
  , "     4.00000     5.00000     6.00000\n"
  , "     7.00000     8.00000     9.00000\n"
  , "    10.00000    11.00000    12.00000\n"
  , "    13.00000    14.00000    15.00000\n"
  , "    16.00000    17.00000    18.00000\n"
  , "    19.00000    20.00000    21.00000\n"
  , "  Hessian / a.u.:\n"
  , "    (6><6)\n"
  , "     [  1.0, 21.0, 21.0, 21.0, 21.0, 21.0\n"
  , "     , 21.0,  2.0, 21.0, 21.0, 21.0, 21.0\n"
  , "     , 21.0, 21.0,  3.0, 21.0, 21.0, 21.0\n"
  , "     , 21.0, 21.0, 21.0,  4.0, 21.0, 21.0\n"
  , "     , 21.0, 21.0, 21.0, 21.0,  5.0, 21.0\n"
  , "     , 21.0, 21.0, 21.0, 21.0, 21.0,  6.0 ]\n"
  , "\n"
  , "#Spicy-Atoms\n"
  , "    C        C0              6     0.50000        0.0009600000     -0.0012400000      0.2814200000         1         2  \n"
  , "    N        N1              7    -0.50000        0.0009600000     -0.0012400000     -0.8757800000         0  \n"
  , "   Fe       Fe2             26     1.40000        0.0004500000     -0.0006100000      2.3073000000         0         3  \n"
  , "    H        H3    P         1     0.40000       -0.0023700000      0.0030100000      4.8693000000         2  \n"
  , "    O        O4              8    -0.80000        2.0592300000     -2.8303200000      3.2801800000         5  \n"
  , "    H        H5              1          No        2.0965600000     -2.8182400000      2.2907100000         4  \n"
  , "    H        H6    P         1     0.40000        2.7085200000     -2.1375800000      3.5614500000  \n"
  ]
