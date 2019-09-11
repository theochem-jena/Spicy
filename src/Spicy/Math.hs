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
{-# LANGUAGE CPP             #-}
{-# LANGUAGE TemplateHaskell #-}
module Spicy.Math
( (<.>)
, vLength
, vDistance
, vAngle
, vCross
, distMat'
, findBondPairs'
) where
import           Data.Array.Accelerate              (Matrix, Vector)
import qualified Data.Array.Accelerate              as A
import qualified Data.Foldable                      as F
import           Data.Sequence                      (Seq)
import qualified Data.Sequence                      as S
import           Prelude                            hiding (cycle, foldl1,
                                                      foldr1, head, init, last,
                                                      maximum, minimum, tail,
                                                      take, takeWhile, (!!))
import qualified Spicy.Math.Internal                as MI
import qualified Spicy.Molecule.Util                as MU
import qualified Data.IntMap                        as IM
import           Data.IntMap                        (IntMap)
import qualified Data.IntSet                        as IS
import           Data.IntSet                        (IntSet)
import           Spicy.Types
#ifdef CUDA
import           Data.Array.Accelerate.LLVM.PTX
#else
import           Data.Array.Accelerate.LLVM.Native
#endif


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

{-|
__PROOF OF CONCEPT FOR ACCELERATE. NOT TO BE TAKEN AS FINAL FUNCION.
-}
distMat' :: Molecule -> Matrix Double
distMat' mol =
  let coordVec = MI.getCoordinates Serial mol
  in  dM coordVec
  where
    dM = $(runQ MI.distMat)

{-|
Build a boolean bond matrix from a Molecule.
-}
findBondPairs' :: Molecule -> IntMap IntSet
findBondPairs' mol =
  let covMat    = cM $ prepareCovalentRadii mol
      coordVec  = MI.getCoordinates Serial mol
      dMat      = dM coordVec
      elemIdxs  = MI.getElementIdxs mol
  in MU.groupTupleSeq $ MI.bondPairsToSeq $ fBP elemIdxs $ bbM dMat covMat
  where cM  = $(runQ $ MI.covRMat 1.3)
        dM  = $(runQ MI.distMat)
        bbM = $(runQ MI.boolBondMatrix)
        fBP = $(runQ MI.findBondPairs)



prepareCovalentRadii :: Molecule -> A.Vector Double
prepareCovalentRadii mol = covRadVec
  where
    msg       = "boolBondMatrix': Something went wrong in looking up covalent radii"
    covRad    = MI.getCovalentRadii mol :: Either String (Vector Double)
    covRadVec = case covRad of
                    Left  _ -> error msg
                    Right c -> c

{-|

-}
-- getBondsIMIS' :: Molecule -> IntMap IntSet
-- getBondsIMIS' mol =
--   let covR = prepareCovalentRadii mol
--       covMat =
--   where

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
