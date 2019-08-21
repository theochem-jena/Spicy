{-|
Module      : Spicy.Molecule.Util
Description : Utilities to manipulate basic molecular data structures.
Copyright   : Phillip Seeber, 2019
License     : GPL-3
Maintainer  : phillip.seeber@uni-jena.de
Stability   : experimental
Portability : POSIX, Windows

This module provides functions to manipulate basic data structures of 'Molecule's, such as indexing.
-}
module Spicy.Molecule.Util
(
) where
import Spicy.Types
import qualified Data.IntMap as IM
import Data.IntMap (IntMap)
import qualified Data.Vector as VB
import Lens.Micro.Platform

{-|
Calculate new indices for 'Molecule's in the 'Individual' scheme, where each fragment is counted by
itself without awareness for other molecules. This relies on sane bonds in the 'Molecule'. If
somewhere bonds are defined to non existing atoms, these bonds will be lost.
-}
-- "RI" = Replaced Indices
reIndexLocal :: Molecule -> Molecule
reIndexLocal mol =
  let molRI        = reIndFrag mol
      molRISubMols = molRI ^. molecule_SubMol
  in  if molRISubMols == VB.empty
        then molRI & molecule_SubMol .~ (VB.map reIndexLocal molRISubMols)
        else molRI
  where
    -- Reindex a layer of 'Molecule' but ignore other layers.
    reIndFrag :: Molecule -> Molecule
    reIndFrag m =
      let atoms      = m ^. molecule_Atoms
          nAtoms     = VB.length nAtoms
          indicesOld = VB.map _atom_Index atoms
          indicesNew = VB.generate nAtoms (\i -> i)
          indReplMap = IM.fromList . VB.toList $ VB.zip indicesOld indicesNew
          bonds      = m ^. molecule_Bonds
          atomsRI    = VB.zipWith (\a i -> a & atom_Index .~ i) atoms indicesNew
          --bondsRI    =

      in  m
            & molecule_Atoms   .~ atomsRI
            & molecule_IndType .~ Individual
