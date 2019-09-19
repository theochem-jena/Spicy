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
, findBonds
, findBondsToGraph
, getBondLength
, calcAnglesBetweenAtoms
) where
import qualified Data.Array.Accelerate             as A
import qualified Data.Foldable                     as F
import qualified Data.Functor                      as FN
import           Data.IntMap                       (IntMap)
import qualified Data.IntMap                       as IM
import           Data.IntSet                       (IntSet)
import           Data.Maybe
import           Data.Sequence                     (Seq)
import qualified Data.Sequence                     as S
-- import           Data.Vector.Unboxed               as VB
import           Data.Graph.Types                  as GT
import           Data.Graph.UGraph                 as UG
import           Prelude                           as P hiding (cycle, foldl1,
                                                    foldr1, head, init, last,
                                                    maximum, minimum, tail,
                                                    take, takeWhile, (!!))
import qualified Spicy.Math.Internal               as MI
import qualified Spicy.Molecule.Util               as MU

import           Spicy.Types
import           Lens.Micro.Platform               as L
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
vAngle a b = acos $ (a <.> b) / (vLength a * vLength b)


{-
Radian to angle conversion on a Seq of angles in rad
-}
radToDegree :: (Floating a) => Seq a -> Seq a
radToDegree = FN.fmap (* ( 180.0 / P.pi ))


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
Quick checking function for testing purposes. Takes a molecule and two indices 
and gives "Maybe" the bond length
-}
getBondLength :: Molecule -> Int -> Int -> Maybe Double
getBondLength mol idx1 idx2 = 
  do
    let atomIM = _molecule_Atoms mol                            :: IntMap Atom
    atom1   <- (^. atom_Coordinates) <$> IM.lookup idx1 atomIM  :: Maybe (Seq Double)
    atom2   <- (^. atom_Coordinates) <$> IM.lookup idx2 atomIM  :: Maybe (Seq Double)

    let bondLength = vDistance atom1 atom2

    return bondLength


{-|
Calculate angles between bonded atoms in a Molecule and give a vector of a triple of atom indices 
and the angle between them.
-}
calcAnglesBetweenAtoms :: Molecule -> UG.UGraph Int () -> Seq (Double, (Int, Int, Int))
calcAnglesBetweenAtoms mol bondGraph = 
  let 
      -- Build primitive Seq's of triples using the bond graph
      tripleSeq  = F.foldl (\ acc (Edge o t _) ->
        let tts       = S.fromList $ GT.adjacentVertices bondGraph t  :: Seq Int
            oos       = S.fromList $ GT.adjacentVertices bondGraph o  :: Seq Int
            trips     = FN.fmap (\tt -> (o, t, tt)) tts               :: Seq (Int, Int, Int)
            altTrips  = FN.fmap (\oo -> (t, o, oo)) oos               :: Seq (Int, Int, Int)
        in  (acc S.><) 
          $ S.filter (\(a, b, c) -> a P./= c P.&& b P./= c) (altTrips S.>< trips)
        ) S.empty $ UG.edges bondGraph
      -- As some triples are duplicate, remove them #Housekeeping :)
      cleanSeq = MU.nubBy (P.==) 
                  $ FN.fmap (\(a, b, c) -> if a P.< c then (a,b,c) else (c,b,a)) tripleSeq
      
      -- Calculate the inner angle of the atom triples
      angles = F.foldl (\acc (a1, a2, a3) -> 
        let atoms     = mol ^. molecule_Atoms                                         :: IntMap Atom
            vec21     = S.zipWith (-) (unsafeCoords atoms a2) (unsafeCoords atoms a1) :: Seq Double 
            vec23     = S.zipWith (-) (unsafeCoords atoms a2) (unsafeCoords atoms a3) :: Seq Double
        in acc S.|> (vAngle vec21 vec23)) S.empty cleanSeq                            :: Seq Double
  in  
      S.zip (radToDegree angles) cleanSeq


{-|
Unsafe retrieving of coordinates. I used the unsafe lookup here, as there is no chance of failure
in the atom indexing (indices are originally drawn from the IntMap, therefore they have to be
"`elem` IntMap.keys")
-}
unsafeCoords :: IntMap Atom -> Int -> Seq Double
unsafeCoords atoms idx = (atoms IM.! idx) ^. atom_Coordinates


{-|
__PROOF OF CONCEPT FOR ACCELERATE. NOT TO BE TAKEN AS FINAL FUNCION.
-}
distMat' :: Molecule -> A.Matrix Double
distMat' mol =
  let coordVec = MI.getCoordinates Serial mol
  in  dM coordVec
  where
#ifdef DEV
    dM = runN MI.distMat
#else
    dM = $(runQ MI.distMat)
#endif
{-|
Accelerate fueled bond finding wrapper function. Gives an IntMap of IntSets with
IntMap indices as bond origins and the IntSet as bond targets for each molecule
-}
findBonds :: Maybe Double -> Molecule -> IntMap IntSet
findBonds covRScaling mol =
  -- Cartesian coordinate vector (size of 3N) of the molecule
  let
      coordVec :: A.Vector Double
      coordVec = MI.getCoordinates Serial mol
      -- Vector of the covalent radii of the atoms in order of the atom indices
      covRadVec :: A.Vector Double
      covRadVec = MI.prepareCovalentRadii mol
      -- Vector of the atom indices
      indices :: A.Vector Int
      indices = MI.getElementIdxs mol
      -- Scalar scaling factor for bond detection using the boolean bond matrix
      scalFac :: A.Scalar Double
      scalFac = A.fromList A.Z [fromMaybe 1.3 covRScaling]

  in  MU.groupTupleSeq . MI.bondPairsToSeq $ bondPairs scalFac indices coordVec covRadVec

  where
    -- Accelerate function chain from the Internal module
#ifdef DEV
    bondPairs = runN MI.accFindBondsChain
#else
    bondPairs = $(runQ MI.accFindBondsChain)
#endif

findBondsToGraph :: Maybe Double -> Molecule -> UGraph Int ()
findBondsToGraph covRScaling mol =
  -- Cartesian coordinate vector (size of 3N) of the molecule
  let coordVec = MI.getCoordinates Serial mol
      -- Vector of the covalent radii of the
      covRadVec = MI.prepareCovalentRadii mol
      indices = MI.getElementIdxs mol
      scalFac = A.fromList A.Z [fromMaybe 1.3 covRScaling]
  in  MI.bondPairsToGraph mol $ bondPairs scalFac indices coordVec covRadVec
  where
#ifdef DEV
    bondPairs = runN MI.accFindBondsChain
#else
    bondPairs = $(runQ MI.accFindBondsChain)
#endif
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
