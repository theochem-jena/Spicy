{-|
Module      : Spicy.Math.Helper
Description : Basic mathematical operations, auxiliary functions for runQ
Copyright   : Phillip Seeber, 2019
License     : GPL-3
Maintainer  : phillip.seeber@uni-jena.de
Stability   : experimental
Portability : POSIX, Windows

Definitions of auxiliary functions to deal with the GHC bug, making it necessary to have all
'A.runQ' expressions as untyped slices. Functions in this module are not meant for direct use, but
wrapped in the 'Spicy.Math' module.
-}
module Spicy.Math.Helper
( vvP
, vmP
, mvP
, mmP
) where
import Data.Array.Accelerate
-- import qualified Data.Array.Accelerate.LLVM.Native           as A
import Data.Array.Accelerate.Numeric.LinearAlgebra
import Prelude hiding ((<>))

{-|
Vector vector dot product.
-}
vvP :: Acc (Vector Double) -> Acc (Vector Double) -> Acc (Scalar Double)
vvP a b = a <.> b

{-|
Vector matrix product.
-}
vmP :: Acc (Vector Double) -> Acc (Matrix Double) -> Acc (Vector Double)
vmP v m = v <# m

{-|
Matrix vector product.
-}
mvP :: Acc (Matrix Double) -> Acc (Vector Double) -> Acc (Vector Double)
mvP m v = m #> v

{-|
Matrix matrix product.
-}
mmP :: Acc (Matrix Double) -> Acc (Matrix Double) -> Acc (Matrix Double)
mmP a b = a <> b
