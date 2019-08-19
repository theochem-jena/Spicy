module Main where
import qualified Data.Array.Accelerate as A
import Spicy.Math


main :: IO ()
main = do
  putStrLn "Executing a simple Accelerate program."
  let vecA = A.fromList (A.Z A.:. 3) [1,2,3]
      vecB = A.fromList (A.Z A.:. 3) [1,2,3]
  print $ vAngle vecA vecB
