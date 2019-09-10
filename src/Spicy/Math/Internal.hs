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
{-# LANGUAGE TypeOperators #-}
module Spicy.Math.Internal
( getCoordinates
, distMat
, covRMat
, boolBondMatrix
, getCovalentRadii
) where
import            Control.Parallel.Strategies
import            Data.Array.Accelerate                         as A
import            Data.Array.Accelerate.Control.Lens
import            Data.Array.Accelerate.IO.Data.Vector.Storable as AVS
import qualified  Data.Foldable                                 as F
import qualified  Data.IntMap                                   as IM
import qualified  Data.Sequence                                 as S
import qualified  Data.Vector.Storable                          as VS
import            Prelude                                       hiding ((/=))
import            Spicy.Types
import qualified  Spicy.Data                                    as D


{-| 
Helper function to get the first element of a triple. Taken from utility-ht == 0.0.14: 
https://hackage.haskell.org/package/utility-ht-0.0.14/docs/src/Data-Tuple-HT.html#fst3
-}
fst3 :: (a,b,c) -> a
fst3 (x,_,_) = x

{-| 
Get the 'Atom' '_atom_Coordinates' from a 'Molecule' and convert to a plain 'VS.Vector'. This is
therefore basically a concatenation of all cartesian coordinates.
-}
getCoordinates :: Strat -> Molecule -> A.Vector Double
getCoordinates strat mol =
  let atomCoords  =
        case strat of
          Serial   ->
            IM.map _atom_Coordinates (mol ^. molecule_Atoms)
          Parallel ->
            IM.map _atom_Coordinates (mol ^. molecule_Atoms) `using` parTraversable rdeepseq
      plainCoords = IM.foldl' (S.><) S.empty atomCoords
      plainVec    = VS.fromList . F.toList $ plainCoords
      vecLength   = VS.length plainVec
  in  AVS.fromVectors (A.Z A.:. vecLength) plainVec

{-|
Calculate the distance matrix from the plain cartesian coordinate vector in \(R^(3 N)\) with \(N\)
being the number of atoms. This functions cannot check, if the cartesian input vector has actually 3
N elements and will crash if the number of elements is not divisable by 3. The calling routine needs
to check this. The resulting matrix will be the square distance matrix.
-}
distMat :: Acc (A.Vector Double) -> Acc (Matrix Double)
distMat v =
  let (Z :. n3) = unlift . shape $ v :: Z :. Exp Int
      n         = n3 `div` 3 :: Exp Int
      dim       = 3 :: Exp Int
      -- Reshape the 3N cartesion input vector to a 3xN matrix with the number of atoms N on
      -- x-axis and (x_n, y_n, z_n) on the y-axis.
      n3Vec     = reshape (lift $ A.Z :. n :. dim) $ v
      -- The x-Axis is now a repetition of the atoms on the y-Axis (which were previously
      -- the x-axis) and z now stores the 3 compotents of the coordinates.
      xVec      = A.replicate (lift $ Z :. n :. All :. All) n3Vec
      -- Transpose the 3D array to swap x- and y-axis and also have the numbersbuildExamples of the atoms on x
      -- again. Strangely the lenses start counting in reverse index order.
      yVec      = transposeOn _2 _3 xVec
  in  -- Overlay the two 3D arrays. The x-y-plane is a table correlating all atom indices with
      -- each other. The z-axis stores the 3 components of the cartesian coordiantes.
      -- Now the 3D arrays will elementwise be subtracted from each other, all results squared,
      -- the z-components folded and square root will be taken and the distance results.
      A.map sqrt . A.sum . A.map (** 2.0) $ A.zipWith (-) xVec yVec



{-|
Retrieve the atom indices as an IntMap for use with the covalentRadii map
-}
getCovalentRadii ::  Molecule -> Either String (A.Vector Double)
getCovalentRadii mol = covRadii
  where 
    -- Get the element numbers in the given molecule
    elementNums = IM.map (fromEnum . _atom_Element) (mol ^. molecule_Atoms )
    -- Get the number of elements in the IntMap
    nrOfAtoms   = IM.size elementNums
    -- Lookup the respective covalent radii from the Spicy.Data.covalentRadiiIM (IntMap)
    -- and assign them to the correspondingbuildExamples atom index
    covRadiiMB  = traverse (`IM.lookup` D.covalentRadiiIM) elementNums
    -- Test if there is some Nothing value in the IntMap (Maybe Double)
    -- If so, give an error message; Else, lift the IntMap out of the Maybe monad -
    -- --> Either String or (IntMap Double) 
    leftMsg     = "Getting the covalent radii failed. Typo in the elements?"
    covRadii    = 
      case covRadiiMB of 
        Nothing   -> Left leftMsg
        Just mat  -> Right $ AVS.fromVectors (Z :. nrOfAtoms) . VS.fromList $ IM.elems mat 
                      

{-|
Calculate the matrix of the sum of the covalent radii for the detection of bonds
--> N x N matrix, symmetric
-}
covRMat :: Exp Double -> Acc (A.Vector Double) -> Acc (A.Matrix Double)
covRMat rFactor covRadii = 
  -- Get the number of atoms from the dimension of the cR vector 
  let (Z :. nrOfAtoms)    = unlift . shape $ covRadii    :: Z :. Exp Int
      -- Build the matrices by replication of the vector in the y direction
      xCovMat             = A.replicate (lift $ Z :. nrOfAtoms :. All) covRadii
      -- The corresponding "y"-matrix is formed by simple transposition the xy plane
      yCovMat             = A.transpose xCovMat  
  -- Get the result by zipping the "3D stack of 2D matrices" using the summation operator
  -- and multiplication by 1.3 (see Literature for the factor)
  -- https://doi.org/10.1063/1.1515483
  in  A.map (* rFactor) $ A.zipWith (+) xCovMat yCovMat  


{-|
Build the boolean bond matrix from the cov
-}
boolBondMatrix :: Acc (A.Matrix Double) -> Acc (A.Matrix Double) -> Acc (A.Matrix Bool)
boolBondMatrix covRMatrix distMatrix = A.zipWith (A.<=) distMatrix covRMatrix


{-| 
Map the boolean bond matrix to an IntMap to get the connectivity and find fragments in the later
run.
-}
findBondPairs :: Molecule -> Acc (A.Matrix Bool) -> Acc (A.Vector (Int, Int))
findBondPairs molA bbMatrix =
    -- Get the ordered element indices in the given molecule
    let atoms       = molA ^. molecule_Atoms
        -- Get the indices from the IntMap of atoms --> [Int]
        idxList     = VS.fromList (IM.keys atoms :: [Int])
        -- Convert to an Acc Vector --> Acc (A.Vector Int)
        elementIdxs = AVS.fromVectors (Z :. IM.size atoms) idxList :: A.Vector Int

        -- Build the index matrix in x-direction
        xMat :: Acc (A.Matrix Int)
        xMat        = A.replicate (lift $ Z :. IM.size atoms :. All) $ A.use elementIdxs 
        
        -- Replicate to get the y-direction
        yMat :: Acc (A.Matrix Int)
        yMat        = A.transpose xMat  

        -- Filter bond pairs using the boolean bond matrix
        (_, o, t)   = A.unzip3 $ (^. _1) $ A.filter (^. _1) $ A.zip3 bbMatrix xMat yMat 

    in  A.zip o t
