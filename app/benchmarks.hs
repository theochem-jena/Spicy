{-
Benchmarks of Spicy functions
-}
{-# LANGUAGE BangPatterns #-}
import           Control.DeepSeq
import           Control.Monad             (replicateM)
import           Criterion.Main
import           Data.Attoparsec.Text.Lazy (many1, parseOnly)
import           Data.List.Split
import qualified Data.Text                 as T
import           Spicy.MolecularSystem
import           Spicy.MolWriter
import           Spicy.Parser
import           Spicy.Types
import           Text.Printf
import qualified Data.Text.IO as T

main = defaultMain
  [ benchmarkGenerators
  , benchmarkParser
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

-- | Benchmarks for generators
benchmarkGenerators = bgroup
  "Generators"
  [ benchmarkTrajGeneratorXYZ
  ]

benchmarkTrajGeneratorXYZ = bgroup
  "XYZ trajectory generator"
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
----------------------------------------------------------------------------------------------------
-- Benchmarks for the parsers
----------------------------------------------------------------------------------------------------
-- | Benchmark for Parsers
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
