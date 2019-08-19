{-|
Module      : Spicy.Math
Description : Basic mathematical operations
Copyright   : Phillip Seeber, 2019
License     : GPL-3
Maintainer  : phillip.seeber@uni-jena.de
Stability   : experimental
Portability : POSIX, Windows

This module defines basic algebraic operations used throughout the program. Numerical heavy and most
other operations are implemented using Accelerate, to provide parallel operations. Note that all
'A.runQ' provided functions must be typed without typeclasses but by concrete types.

The operations here accept some insecurities (like not checking if both vectors of a dot product
have equal lenght) and trust the caller.

-}
{-# LANGUAGE TemplateHaskell, TypeOperators #-}
module Spicy.Math
( (∩)
, (<.>)
, (#>)
, (<#)
, (<>)
, vLength
, vDistance
, vAngle
, vCross
) where
import qualified Data.Array.Accelerate                       as A
import qualified Data.Array.Accelerate.LLVM.Native           as A
import qualified Data.Array.Accelerate.Numeric.LinearAlgebra as A
import           Data.List
import qualified Spicy.Math.Helper as M
import Prelude hiding ((<>))

{-
(<!!>) :: (A.Shape sh, A.Elt e) => A.Array sh e -> Int -> e
arr <!!> ix =
  let accAtIx = A.runQ $ A.unit $ arr A.!! ix :: A.Scalar e
      atIx = head . A.toList $ accAtIx
  in  atIx
  -}



{-|
Intersection (subset) of two lists a and b.
-}
(∩) :: Eq a => [a] -> [a] -> [a]
a ∩ b = a `intersect` b

{-|
Dot product of two 'A.Vector's.
-}
(<.>) :: A.Vector Double -> A.Vector Double -> Double
a <.> b = head . A.toList $ vvP a b
  where
    vvP = $( A.runQ M.vvP )

{-|
'A.Matrix' 'A.Vector' product.
-}
(#>) :: A.Matrix Double -> A.Vector Double -> A.Vector Double
m #> v = mvP m v
  where
    mvP = $( A.runQ M.mvP )

{-|
'A.Vector' 'A.Matrix' product.
-}
(<#) :: A.Vector Double -> A.Matrix Double -> A.Vector Double
v <# m = vmP v m
  where
    vmP = $( A.runQ M.vmP )

{-|
'A.Matrix' 'A.Matrix' product.
-}
(<>) :: A.Matrix Double -> A.Matrix Double -> A.Matrix Double
a <> b = mmP a b
  where
    mmP = $( A.runQ M.mmP )

{-|
Length of a 'A.Vector'.
-}
vLength :: A.Vector Double -> Double
vLength a = head . A.toList $ vL a
  where
    vL = $( A.runQ M.vLength )

{-|
Distance between 2 points ('A.Vector's).
-}
vDistance :: A.Vector Double -> A.Vector Double -> Double
vDistance a b = head . A.toList $ vD a b
  where
    vD = $( A.runQ M.vDistance )

{-|
Angle in radian between 2 'A.Vector's.
-}
vAngle :: A.Vector Double -> A.Vector Double -> Double
vAngle a b = head . A.toList $ vA a b
  where
    vA = $( A.runQ M.vAngle )

{-|
3D cross product of 2 'A.Vectors'
-}
vCross :: A.Vector Double -> A.Vector Double -> A.Vector Double
vCross a b = vC' a b
  where
    vC' = $( A.runQ M.vCross )


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

-- vectorProduct :: A.Vector Double -> A.Vector Double -> A.Scalar Double
-- vectorProduct = $( A.runQ vectorProduct' )
-}

----------------------------------------------------------------------------------------------------
-- Helper functions, not to be exported
{-|
Dimension of a 'A.Vector'.
-}
vDim :: A.Elt a => A.Vector a -> Int
vDim a =
  let A.Z A.:. dim = A.arrayShape a
  in  dim
