{-
Unit testing and golden testing (unit testing based on files) for Spicy. Test all the difficult
parts of Spicy. This is especially Spicy.MolecularSystem and Spicy.Parser.
All tests are required to pass. There is no gray zone!!
-}
{-# LANGUAGE OverloadedStrings #-}
import qualified Data.Array.Accelerate     as A
import           Data.Attoparsec.Text.Lazy
import qualified Data.ByteString.Lazy.UTF8 as LBS
import           Data.Maybe
import qualified Data.Text.Lazy.IO         as T
import qualified Data.Vector               as VB
import qualified Data.Vector.Storable      as VS
import           Spicy.Math
import           Spicy.MolecularSystem
import           Spicy.MolWriter
import           Spicy.Parser
import           Spicy.Types
import           Spicy.UnitTests.Data
import           Test.Tasty
import           Test.Tasty.Golden
import           Test.Tasty.HUnit


-- instance Show Molecule where
--   show = writeSpicy

main :: IO ()
main = defaultMain tests

tests :: TestTree
tests = testGroup "All tests"
  [ testMath
  -- , testParser
  -- , testMolecularSystem
  ]

----------------------------------------------------------------------------------------------------
-- Test cases for math routines
{-|
These test are unit tests for Math parallel and serial computations. Parallel accelerate
calculations are embedded by TemplateHaskell as compiled LLVM code by means of 'A.runQ'. This serves
also as a test for a correct build, as the system is somewhat sensitive with all its dependencies.
-}
testMath :: TestTree
testMath = testGroup "Math"
  [ testDotProduct
  , testVLength
  , testVDistance
  , testVAngle
  , testVCross
  ]

{-|
Dot product of 2 simple 'A.Vector's by Accelerate.
-}
testDotProduct :: TestTree
testDotProduct =
  let vecA = VS.fromList [1, 2 ,3]  :: VS.Vector Double
      vecB = VS.fromList [-7, 8, 9] :: VS.Vector Double
      dotProduct = 36
  in  testCase "Math Dot Product" $
        vecA <.> vecB @?= dotProduct

{-|
Length of a 'A.Vector' by Accelerate.
-}
testVLength :: TestTree
testVLength =
  let vecA = VS.fromList [2, 4, 4] :: VS.Vector Double
      len  = 6
  in  testCase "Vector Length" $
        vLength vecA @?= len

{-|
Distance between 2 points, represented by 'A.Vector's.
-}
testVDistance :: TestTree
testVDistance =
  let vecA = VS.fromList [0, 0, 10] :: VS.Vector Double
      vecB = VS.fromList [0, 0, 0]  :: VS.Vector Double
      dist = 10
  in  testCase "Vectors Distance" $
        vDistance vecA vecB @?= dist

{-|
Anlge between two 'A.Vector's by Accelerate.
-}
testVAngle :: TestTree
testVAngle =
  let vecA  = VS.fromList [0, 0, 1] :: VS.Vector Double
      vecB  = VS.fromList [0, 1, 0] :: VS.Vector Double
      angle = 0.5 * pi
  in  testCase "Vectors Angle" $
        vAngle vecA vecB @?= angle

{-|
Cross product between two 'A.Vector's by Accelerate.
-}
testVCross :: TestTree
testVCross =
  let vecA = VS.fromList [0, 0, 1]  :: VS.Vector Double
      vecB = VS.fromList [0, 1, 0]  :: VS.Vector Double
      vecC = VS.fromList [-1, 0, 0] :: VS.Vector Double
  in  testCase "Vectors Cross Product" $
        vCross vecA vecB @?= vecC

{-
----------------------------------------------------------------------------------------------------
-- Test cases for Parser

{-|
These tests are HUnit tests within the Tasty framework. Correct results and answers are stored in
Spicy.UnitTests.Data, to make sure the parsers work absolutely indepent from environmet. If one of
these tests fail, you are in trouble, as all the others will rely on working parsers and are Golden
Tests instead of UnitTests.
-}
testParser :: TestTree
testParser = testGroup "Parser"
  [ testParserTXYZ1
  , testParserXYZ1
  , testParserMOL21
  , testParserSpicy
  ]

testParserTXYZ1 :: TestTree
testParserTXYZ1 = testCase "Tinker XYZ (1)" $
  (maybeResult $ parse parseTXYZ textHFeCNxH2OTXYZ) @?= Just moleculeHFeCNxH2OTXYZ

testParserXYZ1 :: TestTree
testParserXYZ1 = testCase "Molden XYZ (1)" $
  (maybeResult $ parse parseXYZ testHFeCNxH2OXYZ) @?= Just moleculeHFeCNxH2OXYZ

testParserMOL21 :: TestTree
testParserMOL21 = testCase "SyByl MOL2 (1)" $
 (maybeResult $ parse parseMOL2 testHFeCNxH2OMOL2) @?= Just moleculeHFeCNxH2OMOL2

testParserSpicy :: TestTree
testParserSpicy = testCase "Spicy format (1)" $
  (maybeResult $ parse parseSpicy testHFeCNxH2OSpicy) @?= Just moleculeHFeCNxH2O


----------------------------------------------------------------------------------------------------
-- Test cases for MolecularSystem

testMolecularSystem :: TestTree
testMolecularSystem = testGroup "Molecular System"
  [ testGuessBonds1
  , testGuessBonds2
  , testGuessBonds3
  , testGuessBonds4
  , testGuessBonds5
  , testGuessBonds6
  , testIsolateLayer1
  , testIsolateLayer2
  , testFragmentDetection1
  , testFragmentDetection2
  , testFragmentDetection3
  , testWrapFragmentsToBox1
  , testReplicateSystemAlongAxis1
  , testReplicateSystemAlongAxis2
  , testReplicateSystemAlongAxis3
  , testFindNearestAtom1
  , testFilterByCriteria1
  , testFilterByCriteria2
  , testFilterByCriteria3
  ]

-- Guessing of bonds by colvant radii
testGuessBonds1 :: TestTree
testGuessBonds1 = goldenVsString
  "Guess bonds - defaults (N2 in binding distance)"
  "goldentests/output/N2_bonded__GuessBonds1.txyz" $ do
    raw <- T.readFile "goldentests/input/N2_bonded.xyz"
    case (eitherResult $ parse parseXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let molResult = guessBonds Nothing molInput
        return . LBS.fromString . writeTXYZ $ molResult

testGuessBonds2 :: TestTree
testGuessBonds2 = goldenVsString
  "Guess bonds - defaults (N2 in non-binding distance)"
  "goldentests/output/N2_nonbonded__GuessBonds2.txyz" $ do
    raw <- T.readFile "goldentests/input/N2_nonbonded.xyz"
    case (eitherResult $ parse parseXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let molResult = guessBonds Nothing molInput
        return . LBS.fromString . writeTXYZ $ molResult

testGuessBonds3 :: TestTree
testGuessBonds3 = goldenVsString
  "Guess bonds - custom cutoff (N2 in binding distance)"
  "goldentests/output/N2_bonded__GuessBonds3.txyz" $ do
    raw <- T.readFile "goldentests/input/N2_bonded.xyz"
    case (eitherResult $ parse parseXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let molResult = guessBonds (Just 0.1) molInput
        return . LBS.fromString . writeTXYZ $ molResult

testGuessBonds4 :: TestTree
testGuessBonds4 = goldenVsString
  "Guess bonds - custom cutoff (N2 in non-binding distance)"
  "goldentests/output/N2_nonbonded__GuessBonds4.txyz" $ do
    raw <- T.readFile "goldentests/input/N2_nonbonded.xyz"
    case (eitherResult $ parse parseXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let molResult = guessBonds (Just 7.042254) molInput
        return . LBS.fromString . writeTXYZ $ molResult

testGuessBonds5 :: TestTree
testGuessBonds5 = goldenVsString
  "Guess bonds - 1.2 x R_covalent (Heme like system)"
  "goldentests/output/FePorphyrine__GuessBonds5.txyz" $ do
    raw <- T.readFile "goldentests/input/FePorphyrine.xyz"
    case (eitherResult $ parse parseXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let molResult = guessBonds (Just 1.2) molInput
        return . LBS.fromString . writeTXYZ $ molResult

testGuessBonds6 :: TestTree
testGuessBonds6 = goldenVsString
  "Guess bonds - defaults (sulfate in mixture of H20 and NH3)"
  "goldentests/output/SulfateInSolution__GuessBonds6.txyz" $ do
    raw <- T.readFile "goldentests/input/SulfateInSolution.xyz"
    case (eitherResult $ parse parseXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let molResult = guessBonds Nothing molInput
        return . LBS.fromString . writeTXYZ $ molResult

-- Isolating ONIOM layers and capping dangling bonds
testIsolateLayer1 :: TestTree
testIsolateLayer1 = goldenVsString
  "Isolate ONIOM layer - defaults (Heme like system, isolate Fe-porphyrine)"
  "goldentests/output/FePorphyrine__IsolateLayer1.txyz" $ do
    raw <- T.readFile "goldentests/input/FePorphyrine.txyz"
    case (eitherResult $ parse parseTXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let molResult = isolateLayer [0 .. 35] Nothing Nothing molInput
        return . LBS.fromString . writeTXYZ . fromMaybe moleculeEmpty $ molResult

testIsolateLayer2 :: TestTree
testIsolateLayer2 = goldenVsString
  "Isolate ONIOM layer - fluorine capping, short dinstance (Ru complex with H20 solvent molecules)"
  "goldentests/output/RuKomplex__IsolateLayer2.txyz" $ do
    raw <- T.readFile "goldentests/input/RuKomplex.txyz"
    case (eitherResult $ parse parseTXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let molResult = isolateLayer
              (  [0 .. 19]                                                                -- bPy
              ++ [20]                                                                     -- Ru
              ++ [26, 22, 21, 25, 27, 28, 33, 23, 24, 30, 31, 32, 36, 37, 34, 35, 38, 39] -- bPy
              ++ [46, 42, 44, 41, 43, 45, 49, 47, 48, 50, 51, 52, 54, 55, 57, 56, 58, 53] -- bPy
              ++ [95, 98]                                                                 -- PO
              ++ [126, 127, 128]                                                          -- H2O
              ) (Just F) (Just 0.6) molInput
        return . LBS.fromString . writeTXYZ . fromMaybe moleculeEmpty $ molResult

-- Fragment detection
testFragmentDetection1 :: TestTree
testFragmentDetection1 = goldenVsString
  "Detect fragments - remove all bonds (sulfate in mixture of H20 and NH3)"
  "goldentests/output/SulfateInSolution__FragmentDetection1.xyz" $ do
    raw <- T.readFile "goldentests/input/SulfateInSolution.txyz"
    case (eitherResult $ parse parseTXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let fragments =
              LBS.fromString . concat . map writeTXYZ . snd <$>
              fragmentMolecule RemoveAll molInput
        case fragments of
          Nothing -> return $ LBS.fromString "Failed"
          Just s  -> return s

testFragmentDetection2 :: TestTree
testFragmentDetection2 = goldenVsString
  "Detect fragments - new bond guess (sulfate in mixture of H20 and NH3)"
  "goldentests/output/SulfateInSolution__FragmentDetection2.xyz" $ do
    raw <- T.readFile "goldentests/input/SulfateInSolution.txyz"
    case (eitherResult $ parse parseTXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let fragments =
              LBS.fromString . concat . map writeTXYZ . snd <$>
              fragmentMolecule (NewGuess (Just 1.4)) molInput
        case fragments of
          Nothing -> return $ LBS.fromString "Failed"
          Just s  -> return s

testFragmentDetection3 :: TestTree
testFragmentDetection3 = goldenVsString
  "Detect fragments - keeping supermol bonds (toluene Cl2 mixture periodic)"
  "goldentests/output/TolueneCl2__FragmentDetection3.txyz" $ do
    raw <- T.readFile "goldentests/input/TolueneCl2.txyz"
    case (eitherResult $ parse parseTXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let fragments =
              LBS.fromString . concat . map writeTXYZ . snd <$>
              fragmentMolecule KeepBonds molInput
        case fragments of
          Nothing -> return $ LBS.fromString "Failed"
          Just s  -> return s

-- Wrapping of fragments to unit cell molecule-wise
testWrapFragmentsToBox1 :: TestTree
testWrapFragmentsToBox1 = goldenVsString
  "Wrap molecules to unit cell molecule-wise (toluene Cl2 mixture periodic)"
  "goldentests/output/TolueneCl2__WrapFragmentsToBox1.txyz" $ do
    raw <- T.readFile "goldentests/input/TolueneCl2.txyz"
    case (eitherResult $ parse parseTXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let superMolFragmented =
              LBS.fromString . writeTXYZ . fst . wrapFragmentsToBox (20.0, 20.0, 20.0) <$>
              fragmentMolecule KeepBonds molInput
        case superMolFragmented of
          Nothing -> return $ LBS.fromString "Failed"
          Just s  -> return s

-- | Supercell generation
testReplicateSystemAlongAxis1 :: TestTree
testReplicateSystemAlongAxis1 = goldenVsString
  "Replicate unit cell - x axis (toluene Cl2 mixture periodic)"
  "goldentests/output/TolueneCl2__ReplicateSystemAlongAxis1.xyz" $ do
    raw <- T.readFile "goldentests/input/TolueneCl2.xyz"
    case (eitherResult $ parse parseXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let superCell = replicateSystemAlongAxis (20.0, 20.0, 20.0) AxisX molInput
        return . LBS.fromString . writeXYZ $ superCell

testReplicateSystemAlongAxis2 :: TestTree
testReplicateSystemAlongAxis2 = goldenVsString
  "Replicate unit cell - y axis (toluene Cl2 mixture periodic)"
  "goldentests/output/TolueneCl2__ReplicateSystemAlongAxis2.xyz" $ do
    raw <- T.readFile "goldentests/input/TolueneCl2.xyz"
    case (eitherResult $ parse parseXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let superCell = replicateSystemAlongAxis (20.0, 20.0, 20.0) AxisY molInput
        return . LBS.fromString . writeXYZ $ superCell

testReplicateSystemAlongAxis3 :: TestTree
testReplicateSystemAlongAxis3 = goldenVsString
  "Replicate unit cell - z axis (toluene Cl2 mixture periodic)"
  "goldentests/output/TolueneCl2__ReplicateSystemAlongAxis3.xyz" $ do
    raw <- T.readFile "goldentests/input/TolueneCl2.xyz"
    case (eitherResult $ parse parseXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let superCell = replicateSystemAlongAxis (20.0, 20.0, 20.0) AxisZ molInput
        return . LBS.fromString . writeXYZ $ superCell

-- Nearest neighbour search
testFindNearestAtom1 :: TestTree
testFindNearestAtom1 = goldenVsString
  "Find nearest atom (toluene Cl2 mixture periodic)"
  "goldentests/output/N2_bonded__FindNearestAtom1.dat" $ do
    raw <- T.readFile "goldentests/input/N2_bonded.xyz"
    case (eitherResult $ parse parseXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let nearestInfo = findNearestAtom (0.0, 0.0, 0.5499) molInput
        return . LBS.fromString . show $ nearestInfo

-- Test criterion filtering
testFilterByCriteria1 :: TestTree
testFilterByCriteria1 = goldenVsString
  "Trajectory filtering - distance criterion (azine and phophinin)"
  "goldentests/output/HeteroTraj__FilterByCriteria1.xyz" $ do
    raw <- T.readFile "goldentests/input/HeteroTraj.xyz"
    case (eitherResult $ parse (many1 parseXYZ) raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right trajInput -> do
        let filteredTraj =
              filterByCriteria
              [ fromMaybe False <$> (criterionDistance (2, 16) (< 5.0))
              ] trajInput
        return . LBS.fromString . concat . map writeXYZ $ filteredTraj

testFilterByCriteria2 :: TestTree
testFilterByCriteria2 = goldenVsString
  "Trajectory filtering - 4 atoms angle criterion (azine and phophinin)"
  "goldentests/output/HeteroTraj__FilterByCriteria2.xyz" $ do
    raw <- T.readFile "goldentests/input/HeteroTraj.xyz"
    case (eitherResult $ parse (many1 parseXYZ) raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right trajInput -> do
        let filteredTraj =
              filterByCriteria
              [ fromMaybe False <$> (criterionAngle4Atoms ((5, 2), (2, 16)) (< 1.5708))
              ] trajInput
        return . LBS.fromString . concat . map writeXYZ $ filteredTraj

testFilterByCriteria3 :: TestTree
testFilterByCriteria3 = goldenVsString
  "Trajectory filtering - 2 distance criteria (azine and phophinin)"
  "goldentests/output/HeteroTraj__FilterByCriteria3.xyz" $ do
    raw <- T.readFile "goldentests/input/HeteroTraj.xyz"
    case (eitherResult $ parse (many1 parseXYZ) raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right trajInput -> do
        let filteredTraj =
              filterByCriteria
              [ fromMaybe False <$> (criterionDistance (2, 16) (< 3.85))
              , fromMaybe False <$> (criterionDistance (14, 5) (> 7.85))
              ] trajInput
        return . LBS.fromString . concat . map writeXYZ $ filteredTraj
-}
