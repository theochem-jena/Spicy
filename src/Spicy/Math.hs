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
( (<.>)
, vLength
, vDistance
, vAngle
, vCross
) where
import qualified Data.Foldable as F
import           Data.Sequence (Seq)
import qualified Data.Sequence as S
import           Prelude       hiding (cycle, foldl1, foldr1, head, init, last,
                                maximum, minimum, tail, take, takeWhile, (!!))


{-
(<!!>) :: (VS.Shape sh, VS.Elt e) => VS.Array sh e -> Int -> e
arr <!!> ix =
  let accAtIx = VS.runQ $ VS.unit $ arr VS.!! ix :: VS.Scalar e
      atIx = head . VS.toList $ accAtIx
  in  atIx
  -}

{-|
Dot product of two 'Seq's.
-}
(<.>) :: (Num a) => Seq a -> Seq a -> a
a <.> b = F.sum $ S.zipWith (*) a b

{-|
Length of a 'Seq'.
-}
vLength :: (Floating a) => Seq a -> a
vLength a = sqrt $ a <.> a

{-|
Distance between 2 points ('Seq's).
-}
vDistance :: (Floating a) => Seq a -> Seq a -> a
vDistance a b = vLength $ S.zipWith (-) a b

{-|
Angle in radian between 2 'Seq's.
-}
vAngle :: (Floating a) => Seq a -> Seq a -> a
vAngle a b = acos $ (a <.> b) / ((vLength a) * (vLength b))

{-|
3D cross product of 2 'Seq's.
-}
vCross :: Seq Double -> Seq Double -> Either String (Seq Double)
vCross a b = do
  a1 <- maybeToEither err $ a S.!? 0
  a2 <- maybeToEither err $ a S.!? 1
  a3 <- maybeToEither err $ a S.!? 2
  b1 <- maybeToEither err $ b S.!? 0
  b2 <- maybeToEither err $ b S.!? 1
  b3 <- maybeToEither err $ b S.!? 2
  let c1 = a2 * b3 - a3 * b2
      c2 = a3 * b1 - a1 * b3
      c3 = a1 * b2 - a2 * b1
  return $ S.fromList [c1, c2, c3]
  where
    err = "vCross: Could not get an element from input sequence"

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

-- vectorProduct :: Seq Double -> Seq Double -> VS.Scalar Double
-- vectorProduct = $( VS.runQ vectorProduct' )
-}

----------------------------------------------------------------------------------------------------
-- Helper functions, not to be exported
{-|
Convert a 'Maybe' value to an 'Either' value.
-}
maybeToEither ::
     a          -- ^ 'Left' a will be returned if 'Maybe' was 'Nothing'.
  -> Maybe b    -- ^ 'Right' b will be returned if 'Maybe' was 'Just' b.
  -> Either a b
maybeToEither e Nothing  = Left e
maybeToEither _ (Just a) = Right a
