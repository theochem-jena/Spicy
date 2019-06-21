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

{-
This module defines basic algebraic operations used throughout the program. It
also handles conversion from "R3Vec" to hMatrix's "Vector" type
-}
module Spicy.Math
( (∩)
, r3Vec2hmVec
, hmVec2r3Vec
, hmVecLength
, hmVecDistance
, hmVecAngle
, r3VecNormalVecOfPlane3Points
, r3VecDihedral
) where
import           Data.List
import           Numeric.LinearAlgebra
import           Spicy.Types

-- | intersection (subset) of two lists a and b
(∩) :: Eq a => [a] -> [a] -> [a]
a ∩ b = a `intersect` b

-- | convert from "R3Vec" to "Vector R"
r3Vec2hmVec :: R3Vec -> Vector R
r3Vec2hmVec (a_x, a_y, a_z) = fromList [a_x, a_y, a_z]

-- | convert from "Vector R" to "R3Vec"
hmVec2r3Vec :: Vector R -> Maybe R3Vec
hmVec2r3Vec a
  | (length . toList $ a) == 3 = Just (a_x, a_y, a_z)
  | otherwise = Nothing
  where
    a_x = a `atIndex` 0
    a_y = a `atIndex` 1
    a_z = a `atIndex` 2

-- | calculate the length of a vector
hmVecLength :: Vector R -> R
hmVecLength = sqrt . sum . map (** 2.0) . toList

-- | Distance between 2 points
hmVecDistance :: (Vector R, Vector R) -> R
hmVecDistance (a, b) = hmVecLength $ b - a

-- | Angle between two vectors
hmVecAngle :: (Vector R, Vector R) -> R
hmVecAngle (a, b) = acos $ (a <.> b) / (hmVecLength a * hmVecLength b)

-- | Defines the normal vector of a plane, defined by 3 points
r3VecNormalVecOfPlane3Points :: (Vector R, Vector R, Vector R) -> Vector R
r3VecNormalVecOfPlane3Points (a, b, c) = (b - a) `cross` (c - a)

-- | Dihedral angle between 4 atoms
r3VecDihedral :: (Vector R, Vector R, Vector R, Vector R) -> R
r3VecDihedral (a, b, c, d) = hmVecAngle (p1Normal, p2Normal)
  where
    p1Normal = r3VecNormalVecOfPlane3Points (a, b, c)
    p2Normal = r3VecNormalVecOfPlane3Points (b, c, d)
