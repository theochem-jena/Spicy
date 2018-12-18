import           Data.Attoparsec.Text.Lazy
import           Data.ByteString.Lazy.UTF8      as LBS
import           Data.Either
import           Data.Either.Unwrap
import           Data.Maybe
import           Data.Ord
import qualified Data.Text                 as T
import qualified Data.Text.IO              as T
import           Spicy.MolecularSystem
import           Spicy.MolWriter
import           Spicy.Parser
import           Spicy.Types
import           Spicy.UnitTests.Data
import           System.FilePath
import           System.IO.Unsafe
import           Test.Tasty
import           Test.Tasty.Golden
import           Test.Tasty.HUnit


--main = defaultMain tests
main = defaultMain tests

tests :: TestTree
tests = testGroup "All tests"
  [ testParser
  , testMolecularSystem
  ]

----------------------------------------------------------------------------------------------------
-- Test cases for Parser
----------------------------------------------------------------------------------------------------
-- | These tests are HUnit tests within the Tasty framework. Correct results and answers are stored
-- | in Spicy.UnitTests.Data, to make sure the parsers work absolutely indepent from environmet.
-- | If one of these tests fail, you are in trouble, as all the others will rely on working parsers
-- | and are Golden Tests instead of UnitTests
testParser :: TestTree
testParser = testGroup "Parser"
  [ testParserTXYZ1
  , testParserXYZ1
  , testParserMOL21
  , testParserSpicy
  ]

testParserTXYZ1 = testCase "Tinker XYZ (1)" $
  parseOnly parseTXYZ textHFeCNxH2OTXYZ @?= Right moleculeHFeCNxH2OTXYZ

testParserXYZ1 = testCase "Molden XYZ (1)" $
  parseOnly parseXYZ testHFeCNxH2OXYZ @?= Right moleculeHFeCNxH2OXYZ

testParserMOL21 = testCase "SyByl MOL2 (1)" $
 parseOnly parseMOL2 testHFeCNxH2OMOL2 @?= Right moleculeHFeCNxH2OMOL2

testParserSpicy = testCase "Spicy format (1)" $
  parseOnly parseSpicy testHFeCNxH2OSpicy @?= Right moleculeHFeCNxH2O

----------------------------------------------------------------------------------------------------
-- Test cases for MolecularSystem
----------------------------------------------------------------------------------------------------
testMolecularSystem :: TestTree
testMolecularSystem = testGroup "Molecular System"
  [ testGuessBonds1
  , testGuessBonds2
  , testGuessBonds3
  , testGuessBonds4
  , testGuessBonds5
  , testGuessBonds6
  ]

-- | Guessing of bonds
testGuessBonds1 = goldenVsString
  "Guess bonds - defaults (N2 in binding distance)"
  "goldentests/output/N2_bonded__GuessBonds1.txyz" $ do
    raw <- T.readFile "goldentests/input/N2_bonded.xyz"
    case (parseOnly parseXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let molResult = guessBonds Nothing molInput
        return . LBS.fromString . writeTXYZ $ molResult

testGuessBonds2 = goldenVsString
  "Guess bonds - defaults (N2 in non-binding distance)"
  "goldentests/output/N2_nonbonded__GuessBonds2.txyz" $ do
    raw <- T.readFile "goldentests/input/N2_nonbonded.xyz"
    case (parseOnly parseXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let molResult = guessBonds Nothing molInput
        return . LBS.fromString . writeTXYZ $ molResult

testGuessBonds3 = goldenVsString
  "Guess bonds - custom cutoff (N2 in binding distance)"
  "goldentests/output/N2_bonded__GuessBonds3.txyz" $ do
    raw <- T.readFile "goldentests/input/N2_bonded.xyz"
    case (parseOnly parseXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let molResult = guessBonds (Just 0.1) molInput
        return . LBS.fromString . writeTXYZ $ molResult

testGuessBonds4 = goldenVsString
  "Guess bonds - custom cutoff (N2 in non-binding distance)"
  "goldentests/output/N2_nonbonded__GuessBonds4.txyz" $ do
    raw <- T.readFile "goldentests/input/N2_nonbonded.xyz"
    case (parseOnly parseXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let molResult = guessBonds (Just 7.042254) molInput
        return . LBS.fromString . writeTXYZ $ molResult

testGuessBonds5 = goldenVsString
  "Guess bonds - 1.2 x R_covalent (Heme like system)"
  "goldentests/output/FePorphyrine__GuessBonds5.txyz" $ do
    raw <- T.readFile "goldentests/input/FePorphyrine.xyz"
    case (parseOnly parseXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let molResult = guessBonds (Just 1.2) molInput
        return . LBS.fromString . writeTXYZ $ molResult

testGuessBonds6 = goldenVsString
  "Guess bonds - defaults (sulfate in mixture of H20 and NH3)"
  "goldentests/output/SulfateInSolution__GuessBonds6.txyz" $ do
    raw <- T.readFile "goldentests/input/SulfateInSolution.xyz"
    case (parseOnly parseXYZ raw) of
      Left _ -> return $ LBS.fromString "Failed"
      Right molInput -> do
        let molResult = guessBonds Nothing molInput
        return . LBS.fromString . writeTXYZ $ molResult
