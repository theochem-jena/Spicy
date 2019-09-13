module Main where

import Spicy.Math as SM
-- import Spicy.Types
import Spicy.Parser
-- import Data.Array.Accelerate as A
-- import Data.Array.Accelerate.LLVM.Native as LLVM
-- import Data.Graph.UGraph as UG
import Debug.Trace
import Data.Graph.Visualize as GV
import Data.Text.Lazy.IO as TIO
import Data.Attoparsec.Text.Lazy (parse, Result(Fail, Done))

main :: IO ()
main = do
  fileContent <- TIO.readFile "goldentests/input/mDCB+.xyz"
  info        <- case parse parseXYZ fileContent of
    Fail _ _ e  -> error e
    Done _ mol  -> do
      let graph = SM.findBondsToGraph Nothing mol
      _ <- GV.plotUGraphPng graph "test"
      traceShowM graph

  print info
