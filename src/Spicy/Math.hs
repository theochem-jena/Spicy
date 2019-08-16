{-|
Module      : Spicy.Math
Description : Basic mathematical operations
Copyright   : Phillip Seeber, 2019
License     : GPL-3
Maintainer  : phillip.seeber@uni-jena.de
Stability   : experimental
Portability : POSIX, Windows

This module defines basic algebraic operations used throughout the program. It
also handles conversion from "R3Vec" to hMatrix's "Vector" type

    * highest level region has highest index

    * lowest level region has index 0 and contains the complete system
-}
{-# LANGUAGE TemplateHaskell #-}
module Spicy.Math
( (∩)
) where
import qualified Data.Array.Accelerate                       as A
import qualified Data.Array.Accelerate.LLVM.Native           as A
import qualified Data.Array.Accelerate.Numeric.LinearAlgebra as A
import           Data.List
--import           Spicy.MathHelper
--import           Spicy.Types

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
(<.>) :: A.Numeric a => A.Vector a -> A.Vector a ->  a
a <.> b = head . A.toList $ (A.runN (A.<.>)) a b

{-|
Length of a 'A.Vector'.
-}
vLength :: (A.Numeric a, Floating a) => A.Vector a -> a
vLength a = sqrt (a <.> a)

{-|
Distance between 2 points ('A.Vector's).
-}
-- For Accelerate's fusion, don't call 'vLength' here.
vDistance :: (A.Numeric a, Floating a) => A.Vector a -> A.Vector a -> a
vDistance a b = sqrt . head . A.toList $ (A.runN f) a b
  where
    f :: (A.Numeric a, Floating a) => A.Acc (A.Vector a) -> A.Acc (A.Vector a) -> A.Acc (A.Scalar a)
    f v1 v2 = (\x -> x A.<.> x) $ A.zipWith (-) v1 v2


{-|
Angle between to 'A.Vector's in radian.
-}
-- For Accelerate's fusion, don't call other 'Spicy.Math' functions here.
vAngle :: A.Vector Double -> A.Vector Double -> Double
vAngle a b = head . A.toList $ (A.runN f) a b
  where
    f :: A.Acc (A.Vector Double) -> A.Acc (A.Vector Double) -> A.Acc (A.Scalar Double)
    f x y = A.zipWith (/) (dividend x y) (divisor x y)-- (\x -> x A.<.> x) $ A.zipWith (-) v1 v2
    --
    dividend :: A.Acc (A.Vector Double) -> A.Acc (A.Vector Double) -> A.Acc (A.Scalar Double)
    dividend x y = x A.<.> y
    --
    divisor :: A.Acc (A.Vector Double) -> A.Acc (A.Vector Double) -> A.Acc (A.Scalar Double)
    divisor x y = A.zipWith (*) (vLength' x) (vLength' y)
    --
    vLength' :: A.Acc (A.Vector Double) -> A.Acc (A.Scalar Double)
    vLength' x = A.map A.sqrt $ (x A.<.> x)


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
-}

-- vectorProduct :: A.Vector Double -> A.Vector Double -> A.Scalar Double
-- vectorProduct = $( A.runQ vectorProduct' )
