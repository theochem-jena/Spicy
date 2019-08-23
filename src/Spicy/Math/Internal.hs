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
module Spicy.Math.Internal
( getCoordinates1
, getCoordinates2
) where
import qualified Data.Array.Accelerate                       as A
import           Lens.Micro.Platform
import           Spicy.Types
import qualified Data.Vector           as VB
import qualified Data.Vector.Storable  as VS
import Control.Parallel.Strategies
import qualified Data.IntMap as IM
import Data.IntMap (IntMap)

{-|
Get the 'Atom' '_atom_Coordinates' from a 'Molecule' and convert to a plain 'VS.Vector'. This is therefore
basically a concatenation of all cartesian coordinates
-}
getCoordinates1 :: Strat -> Molecule -> VS.Vector Double
getCoordinates1 strat mol =
  let atomCoords  = case strat of
        Serial   -> IM.map _atom_Coordinates (mol ^. molecule_Atoms)
        Parallel -> IM.map _atom_Coordinates (mol ^. molecule_Atoms) `using` parTraversable rdeepseq
      plainCoords = IM.foldl' (VS.++) VS.empty atomCoords
  in  plainCoords

getCoordinates2 :: Strat -> Molecule -> VS.Vector Double
getCoordinates2 s m = VS.concat $ m ^.. molecule_Atoms . each . atom_Coordinates
