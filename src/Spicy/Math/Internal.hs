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
( 
) where
import           Data.Array.Accelerate                       as A
import           Data.Array.Accelerate.Numeric.LinearAlgebra as A
import           Lens.Micro.Platform
import           Spicy.Types
import qualified Data.Vector           as VB
import qualified Data.Vector.Storable  as VS
