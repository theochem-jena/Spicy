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
, vLength
, vDistance
, vAngle
, vCross
) where
import           Data.Array.Accelerate                       as A
import           Data.Array.Accelerate.Numeric.LinearAlgebra as A


{-|
'A.Vector' 'A.Vector' dot product.
-}
vvP :: A.Acc (A.Vector Double) -> A.Acc (A.Vector Double) -> A.Acc (Scalar Double)
vvP a b = a A.<.> b

{-|
'A.Vector' 'A.Matrix' product.
-}
vmP :: A.Acc (A.Vector Double) -> A.Acc (A.Matrix Double) -> A.Acc (A.Vector Double)
vmP v m = v A.<# m

{-|
'A.Matrix' 'A.Vector' product.
-}
mvP :: A.Acc (A.Matrix Double) -> A.Acc (A.Vector Double) -> A.Acc (A.Vector Double)
mvP m v = m A.#> v

{-|
'A.Matrix' 'A.Matrix' product.
-}
mmP :: A.Acc (A.Matrix Double) -> A.Acc (A.Matrix Double) -> A.Acc (A.Matrix Double)
mmP a b = a A.<> b

{-|
Length of a 'A.Vector'.
-}
vLength :: A.Acc (A.Vector Double) -> A.Acc (A.Scalar Double)
vLength a = A.map sqrt $ a <.> a

{-|
Distance between 2 'A.Vector's.
-}
vDistance :: A.Acc (A.Vector Double) -> A.Acc (A.Vector Double) -> A.Acc (A.Scalar Double)
vDistance a b = vLength $ A.zipWith (-) a b

{-|
Angle in radian between 2 'A.Vector's.
-}
vAngle :: A.Acc (A.Vector Double) -> A.Acc (A.Vector Double) -> A.Acc (A.Scalar Double)
vAngle a b =
  let dividend = a A.<.> b
      divisor  = A.zipWith (*) (vLength a) (vLength b)
  in  A.map A.acos $ A.zipWith (/) dividend divisor

vCross :: A.Acc (A.Vector Double) -> A.Acc (A.Vector Double) -> A.Acc (A.Vector Double)
vCross a b =
  let a1 = a A.!! 0
      a2 = a A.!! 1
      a3 = a A.!! 2
      b1 = b A.!! 0
      b2 = b A.!! 1
      b3 = b A.!! 2
      c1 = a2 * b3 - a3 * b2
      c2 = a3 * b1 - a1 * b3
      c3 = a1 * b2 - a2 * b1
  in  A.generate (A.index1 3) (\ix ->
        let i = A.unindex1 ix
        in  A.caseof i
              [ ((\j -> j A.== 0), c1)
              , ((\j -> j A.== 1), c2)
              , ((\j -> j A.== 2), c3)
              ]
              0
      )
