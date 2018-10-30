{-
This module defines basic algebraic operations used throughout the program. It
also handles conversion from "R3Vec" to hMatrix's "Vector" type
-}
module Spicy.Math
( (∩)
, r3Vec2hmVec
, hmVec2r3Vec
, hmVecLength
) where
import           Data.List
import           Numeric.LinearAlgebra
import           Numeric.LinearAlgebra.Data
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
hmVecLength :: Vector R -> Double
hmVecLength = sqrt . sum . map (^2) . toList
