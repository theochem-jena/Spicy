import           Data.Attoparsec.Text.Lazy
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
import           System.IO.Unsafe
import           Test.Tasty
import           Test.Tasty.HUnit

instance Show Molecule where
  show = write_MOL2

main = defaultMain tests

tests :: TestTree
tests = testGroup "All tests"
  [ testParser
  ]

----------------------------------------------------------------------------------------------------
-- Test cases for Parser
----------------------------------------------------------------------------------------------------
testParser :: TestTree
testParser = testGroup "Parser"
  [ testParserTXYZ1
  , testParserXYZ1
  , testParserMOL21
  ]

testParserTXYZ1 = testCase "Tinker XYZ (1)" $
  parseOnly parse_TXYZ textHFeCNxH2OTXYZ @?= Right moleculeHFeCNxH2OTXYZ

testParserXYZ1 = testCase "Molden XYZ (1)" $
  parseOnly parse_XYZ testHFeCNxH2OXYZ @?= Right moleculeHFeCNxH2OXYZ

testParserMOL21 = testCase "SyByl MOL2 (1)" $
 parseOnly parse_MOL2 testHFeCNxH2OMOL2 @?= Right moleculeHFeCNxH2OMOL2

----------------------------------------------------------------------------------------------------
-- Test cases for MolecularSystem
----------------------------------------------------------------------------------------------------
{-
testMolecularSystem :: TestTree
testMolecularSystem = testGroup "Molecular System"
  [ testMolecularSystemBondGuess1
  ]
testMolecularSystemBondGuess1 = testCase "Guessing bonds" $
  guessBonds (Just 1.4) mol1_xyz @?= moleculeHFeCNxH2O
-}
