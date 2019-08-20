{-|
Module      : Spicy.Math
Description : Basic mathematical operations
Copyright   : Phillip Seeber, 2019
License     : GPL-3
Maintainer  : phillip.seeber@uni-jenVS.de
Stability   : experimental
Portability : POSIX, Windows

This module defines basic algebraic operations used throughout the program. Numerical heavy and most
other operations are implemented using Accelerate, to provide parallel operations. Note that all
'VS.runQ' provided functions must be typed without typeclasses but by concrete types.

The operations here accept some insecurities (like not checking if both vectors of a dot product
have equal lenght) and trust the caller.

-}
module Spicy.Math
( (∩)
, (<.>)
, vLength
, vDistance
, vAngle
, vCross
) where
import qualified Data.Array.Accelerate                       as A
import qualified Data.Array.Accelerate.LLVM.Native           as A
import qualified Data.Array.Accelerate.Numeric.LinearAlgebra as A
import           Data.List
import           Prelude                                     hiding ((<>))
import qualified Spicy.Math.Internal                         as M
import qualified Data.Vector           as VB
import qualified Data.Vector.Storable  as VS

{-
(<!!>) :: (VS.Shape sh, VS.Elt e) => VS.Array sh e -> Int -> e
arr <!!> ix =
  let accAtIx = VS.runQ $ VS.unit $ arr VS.!! ix :: VS.Scalar e
      atIx = head . VS.toList $ accAtIx
  in  atIx
  -}



{-|
Intersection (subset) of two lists a and b.
-}
(∩) :: Eq a => [a] -> [a] -> [a]
a ∩ b = a `intersect` b

{-|
Dot product of two 'VS.Vector's.
-}
(<.>) :: (VS.Storable a, Num a) => VS.Vector a -> VS.Vector a -> a
a <.> b = VS.sum $ VS.zipWith (*) a b

{-|
Length of a 'VS.Vector'.
-}
vLength :: (VS.Storable a, Floating a) => VS.Vector a -> a
vLength a = sqrt $ a <.> a

{-|
Distance between 2 points ('VS.Vector's).
-}
vDistance :: (VS.Storable a, Floating a) => VS.Vector a -> VS.Vector a -> a
vDistance a b = vLength $ VS.zipWith (-) a b

{-|
Angle in radian between 2 'VS.Vector's.
-}
vAngle :: (VS.Storable a, Floating a) => VS.Vector a -> VS.Vector a -> a
vAngle a b = acos $ (a <.> b) / ((vLength a) * (vLength b))

{-|
3D cross product of 2 'VS.Vectors'
-}
vCross :: VS.Vector Double -> VS.Vector Double -> VS.Vector Double
vCross a b =
  let a1 = a VS.! 0
      a2 = a VS.! 1
      a3 = a VS.! 2
      b1 = b VS.! 0
      b2 = b VS.! 1
      b3 = b VS.! 2
      c1 = a2 * b3 - a3 * b2
      c2 = a3 * b1 - a1 * b3
      c3 = a1 * b2 - a2 * b1
  in  VS.fromList [c1, c2, c3]




{-|

-}


{-
-- | Defines the normal vector of a plane, defined by 3 points
r3VecNormalVecOfPlane3Points :: (Vector R, Vector R, Vector R) -> Vector R
r3VecNormalVecOfPlane3Points (a, b, c) = (b - a) `cross` (c - a)

-- | Dihedral angle between 4 atoms
r3VecDihedral :: (Vector R, Vector R, Vector R, Vector R) -> R
r3VecDihedral (a, b, c, d) = hmVecAngle (p1Normal, p2Normal)
  where
    p1Normal = r3VecNormalVecOfPlane3Points (a, b, c)
    p2Normal = r3VecNormalVecOfPlane3Points (b, c, d)

-- vectorProduct :: VS.Vector Double -> VS.Vector Double -> VS.Scalar Double
-- vectorProduct = $( VS.runQ vectorProduct' )
-}

----------------------------------------------------------------------------------------------------
-- Helper functions, not to be exported

{-|
Convert
-}

{-
parVector :: NFData a => Strategy (V.Vector a)
parVector vec =
  let vLen = V.length vec
      half = vLen `div` 2
      minChunk = 10
  in  if vLen > minChunk
      then do
        let v1 = V.unsafeSlice 0 half vec
            v2 = V.unsafeSlice half (vLen - half) vec
        parVector v1
        parVector v2
        return vec
      else
        evalChunk (vLen-1) >>
        return vec
  where
  evalChunk 0 = rpar (rdeepseq (vec V.! 0)) >> return vec
  evalChunk i = rpar (rdeepseq (vec V.! i)) >> evalChunk (i-1)
-}
