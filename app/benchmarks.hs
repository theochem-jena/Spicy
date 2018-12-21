{-
Benchmarks of Spicy functions
-}
{-# LANGUAGE BangPatterns #-}
import           Control.Monad             (replicateM)
import           Criterion.Main
import           Data.Attoparsec.Text.Lazy (many1, parseOnly)
import           Data.IntSet               (IntSet)
import qualified Data.IntSet               as I
import           Data.List.Split
import qualified Data.Text                 as T
import qualified Data.Text.IO              as T
import           Lens.Micro.Platform
import           Spicy.MolecularSystem
import           Spicy.MolWriter
import           Spicy.Parser
import           Spicy.Types
import           Text.Printf

main = defaultMain
  [ benchmarkGenerators
  , benchmarkParser
  , benchmarkMolecularSystem
  ]
----------------------------------------------------------------------------------------------------
-- Generate test data
----------------------------------------------------------------------------------------------------
-- | generate a very simple repetitive trajectory in XYZ format
generateTrajectoryXYZ :: Int -> Int -> T.Text
generateTrajectoryXYZ nAtoms nFrames = T.pack frames
  where
    atomLine = "Rb   1.00000  -15.031654 0.00354\n"
    frame =
      (show nAtoms) ++ "\n" ++
      "Comment\n" ++
      ( concat $
        replicate nAtoms atomLine
      )
    frames =
      concat $
      replicate nFrames frame

generateNonbondedMolecules :: Int -> Molecule
generateNonbondedMolecules nMols = superMolecule
  where
    molecule = Molecule
      { _molecule_Label = "Phosphinin"
      , _molecule_Atoms = atomsCarbons ++ atomsHydrogen ++ atomsPhosphor
      , _molecule_Energy = Nothing
      , _molecule_Gradient = Nothing
      , _molecule_Hessian = Nothing
      }
    templateAtom = Atom
        { _atom_Element = undefined
        , _atom_Label = ""
        , _atom_IsPseudo = False
        , _atom_FFType = ""
        , _atom_PCharge = Nothing
        , _atom_Coordinates = undefined
        , _atom_Connectivity = I.empty
        }
    atomsCarbons = zipWith (\a c -> a & atom_Coordinates .~ c & atom_Element .~ C)
      (replicate 5 templateAtom)
      [ ( 1.18770,  0.66347,  2.50000)
      , ( 1.20847,  2.06088,  2.50000)
      , (-1.18770,  0.66348,  2.50000)
      , ( 0.00000,  2.76139,  2.50000)
      , (-1.20847,  2.06088,  2.50000)
      ]
    atomsHydrogen = zipWith (\a c -> a & atom_Coordinates .~ c & atom_Element .~ H)
      (replicate 5 templateAtom)
      [ (-2.15024,  2.59423,  2.50000)
      , ( 0.00000,  3.84372,  2.50000)
      , (-2.11521,  0.10669,  2.50000)
      , ( 2.11521,  0.10669,  2.50000)
      , ( 2.15024,  2.59423,  2.50000)
      ]
    atomsPhosphor = zipWith (\a c -> a & atom_Coordinates .~ c & atom_Element .~ P)
      (replicate 1 templateAtom)
      [ (0.00000,  0.00000,  2.50000)
      ]
    shiftVec =
      [ (\(x, y, z) -> (6.0 * fromIntegral x, 6.0 * fromIntegral y, 3.0 *fromIntegral z)) i
      | i <- (\x y z -> (x, y, z)) <$>
          [1 .. l] <*>
          [1 .. l] <*>
          [1 .. l]
      ]
      where
        l = ceiling . (**(1/3)) . fromIntegral $ nMols
    superMolecule =
      molecule & molecule_Atoms .~
      ( concat . map _molecule_Atoms
      $ zipWith (shiftFragment) (take nMols shiftVec) (repeat molecule)
      )
--shiftVec :: Int -> [(Double, Double, Double)]
shiftVec n =
  take n $
  [ (\(x, y, z) -> (x, y, z)) i
  | i <- (\x y z -> (x, y, z)) <$>
      [1 .. l] <*>
      [1 .. l] <*>
      [1 .. l]
  ]
  where
    l = ceiling . (**(1/3)) . fromIntegral $ n

-- | Benchmarks for generators
benchmarkGenerators = bgroup
  "Generators"
  [ benchmarkTrajGeneratorXYZ
  , benchmarkGenerateNonbondedMolecules
  ]

benchmarkTrajGeneratorXYZ = bgroup
  "XYZ trajectory"
  [ benchmarkTrajGeneratorXYZ100A100F
  , benchmarkTrajGeneratorXYZ100A500F
  , benchmarkTrajGeneratorXYZ100A1000F
  , benchmarkTrajGeneratorXYZ100A5000F
  ]

benchmarkTrajGeneratorXYZ100A100F = bench
  "100 atoms, 100 frames" $
  nf (generateTrajectoryXYZ 100) 100

benchmarkTrajGeneratorXYZ100A500F = bench
  "100 atoms, 500 frames" $
  nf (generateTrajectoryXYZ 100) 500

benchmarkTrajGeneratorXYZ100A1000F = bench
  "100 atoms, 1000 frames" $
  nf (generateTrajectoryXYZ 100) 1000

benchmarkTrajGeneratorXYZ100A5000F = bench
  "100 atoms, 5000 frames" $
  nf (generateTrajectoryXYZ 100) 5000

benchmarkGenerateNonbondedMolecules = bgroup
  "Multiple molecules"
  [ benchmarkGenerateNonbondedMoleculesSmall
  , benchmarkGenerateNonbondedMoleculesMedium
  , benchmarkGenerateNonbondedMoleculesLarge
  ]

benchmarkGenerateNonbondedMoleculesSmall = bench
  "10 molecules" $
  nf generateNonbondedMolecules 10

benchmarkGenerateNonbondedMoleculesMedium = bench
  "50 molecules" $
  nf generateNonbondedMolecules 50

benchmarkGenerateNonbondedMoleculesLarge = bench
  "100 molecules" $
  nf generateNonbondedMolecules 100

----------------------------------------------------------------------------------------------------
-- Benchmarks for the parsers
----------------------------------------------------------------------------------------------------
benchmarkParser = bgroup
  "Parser"
  [ benchmarkParserXYZ
  ]

benchmarkParserXYZ = bgroup
  "XYZ Trajectory"
  [ benchmarkParserXYZ100A100F
  , benchmarkParserXYZ100A1000F
  , benchmarkParserXYZ1000A100F
  ]

benchmarkParserXYZ100A100F = bench
  "100 atoms, 100 frames" $
  nf (parseOnly (many1 parseXYZ)) (generateTrajectoryXYZ 100 100)

benchmarkParserXYZ100A1000F = bench
  "100 atoms, 1000 frames" $
  nf (parseOnly (many1 parseXYZ)) (generateTrajectoryXYZ 100 1000)

benchmarkParserXYZ1000A100F = bench
  "1000 atoms, 100 frames" $
  nf (parseOnly (many1 parseXYZ)) (generateTrajectoryXYZ 1000 100)

-----------------------------------------------------------------------------------------------------
-- Benchmarks MolecularSystem
----------------------------------------------------------------------------------------------------

benchmarkMolecularSystem = bgroup
  "Molecular system"
  [ benchmarkGuessBonds
  ]

benchmarkGuessBonds = bgroup
  "Guess bonds"
  [ benchmarkGuessBondsSmall
  , benchmarkGuessBondsMedium
  , benchmarkGuessBondsLarge
  ]

benchmarkGuessBondsSmall = bench
  "Small (10 molecules)" $
  nf (guessBonds (Just 1.2)) (generateNonbondedMolecules 10)

benchmarkGuessBondsMedium = bench
  "Medium (50 molecules)" $
  nf (guessBonds (Just 1.2)) (generateNonbondedMolecules 50)

benchmarkGuessBondsLarge = bench
  "Medium (100 molecules)" $
  nf (guessBonds (Just 1.2)) (generateNonbondedMolecules 100)
