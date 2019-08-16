module Spicy.MathHelper
( vectorProduct'
) where
import           Data.List
import qualified Data.Array.Accelerate as A
import qualified Data.Array.Accelerate.LLVM.Native as A
import qualified Data.Array.Accelerate.Numeric.LinearAlgebra as A

import           Spicy.Types

vectorProduct' :: A.Acc (A.Vector Double) -> A.Acc (A.Vector Double) -> A.Acc (A.Scalar Double)
vectorProduct' = (A.<.>)
