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
import           Data.IntMap         (IntMap, (!?))
import qualified Data.IntMap         as IM
import           Data.IntSet         (IntSet, (\\))
import qualified Data.IntSet         as IS
import           Data.Maybe
import qualified Data.Vector         as VB
import           Lens.Micro.Platform
import           Spicy.Types

{-
"IS" = IntSet
"IM" = IntMap
"oK" = old key
"nK" = new key
"RI" = replaced index
"oS" = old set
"nS" = new set
"oA" = old atom
"nA" = new atom
-}

{-|
Reindex a complete 'Molecule', including all its deeper layers in '_molecule_SubMol') by mappings
from a global replacement Map, mapping old to new indices. This function assumes, that your molecule
is sane in the overall assumptions of this program. This means that the lower layers obey the
counting scheme of the atoms of the higher layers and pseudoatoms come last.
-}
reIndexMolecule :: IntMap Int -> Molecule -> Either String Molecule
reIndexMolecule repMap mol = do
  -- Get the submolecules
  let subMols = mol ^. molecule_SubMol
  -- Reindex the current layer only
  molRI <- reIndexMoleculeLayer repMap mol
    -- If the molecule has no submolecules:
  if VB.null subMols
    -- Then we are done.
    then Right molRI
    -- Else we need to reindex the deeper layers also.
    else do
      subMolsRI <- VB.sequence . VB.map (reIndexMolecule repMap) $ subMols
      reIndexMolecule repMap $ molRI & molecule_SubMol .~ subMolsRI

{-|
Reindex the '_molecule_Atoms' and '_molecule_Bonds' of a single layer of a molecule (ignoring
anything in the '_molecule_SubMol' field). While the completeness of the reindexing is checked and
incompleteness of the replacement 'IntMap' 'Int' will result in 'Left' 'String', it is not checked
if the 'Atom''s indexing is sane and indices are unique.
-}
reIndexMoleculeLayer ::
     IntMap Int             -- ^ 'IntMap' with mappings from old indices to new indices (bonds and
                            --   atoms).
  -> Molecule               -- ^ Molecule to reindex.
  -> Either String Molecule -- ^ 'Left' 'String' in case the mapping from old keys to new keys is
                            --   incomplete.
reIndexMoleculeLayer repMap mol
  -- If old Atom indices cannot be completely replaced by the Map:
  | not $ checkRepMapCompleteIS repMap (IM.keysSet $ mol ^. molecule_Atoms) =
      Left "reIndexMoleculeLayer: The remapping of indices is not complete for the atom indices."
  -- If old bond type data cannot be completely replaced by the Map:
  | not $ checkRepMapCompleteIMIS repMap (mol ^. molecule_Bonds)            =
      Left "reIndexMoleculeLayer: The remapping of indices is not complete for the bond data."
  -- If both checks went fine:
  | otherwise                                                               =
      Right $ mol
        -- Use the (%~) lens to update the atoms indices of a molecule with new indices.
        & molecule_Atoms %~ replaceIMKeys repMap
        -- Use the (%~) lens to update keys and values (IntSet) of the bond type data structure.
        & molecule_Bonds %~ replaceIMIS repMap

{-|
Check if 'IntMap' is complete to replace all values in an 'IntSet' (any old 'IS.Key' can be replaced
by a new 'IS.Key').  Gives 'True' if complete, 'False' if the replacement 'IntMap' 'Int' has holes.
-}
checkRepMapCompleteIS ::
     IntMap Int -- ^ 'IntMap' containing the mapping from old 'IS.Key's to new 'IS.Key's.
  -> IntSet     -- ^ 'IntSet' in which values should be replaced.
  -> Bool       -- ^ Result.
checkRepMapCompleteIS repMap is =
  let -- Keys, that can be replaced
      repKeys  = IM.keysSet repMap
      -- Keys that shall be replaced minus keys that can be replaced
      lostKeys = is \\ repKeys
  in  IS.null lostKeys

{-|
Check if 'IntMap' is complete to replace all 'IM.Key's from the old 'IntMap' by new 'IM.Key's. Gives
'True' if the 'IntMap' with replacements is complete and 'False' otherwise.
-}
checkRepMapCompleteIM ::
     IntMap Int -- ^ 'IntMap' containing the mapping from old 'IS.Key's to new 'IS.Key's.
  -> IntMap a   -- ^ 'IntMap' in which 'IM.Key's should be replaced.
  -> Bool       -- ^ Result.
checkRepMapCompleteIM repMap im =
  let -- Keys, that can be replaced
      repKeys = IM.keysSet repMap
      -- Keys that shall be replaced
      oldKeys = IM.keysSet im
      -- Keys that cannot be replaced
      lostKeys = oldKeys \\ repKeys
  in  IS.null lostKeys


{-|
Check if 'IntMap' is complete to replace all values in an 'IntMap' 'InSet' type construction
(replacing both the lookup keys in the 'IntMap', as well as all values in the 'IntSet'). Gives
'True' if complete and 'False' otherwise.
-}
checkRepMapCompleteIMIS ::
     IntMap Int    -- ^ 'IntMap' containing the mapping from old 'IS.Key's to new 'IS.Key's.
  -> IntMap IntSet -- ^ 'IntMap' 'IntSet' in which values and keys should be replaced.
  -> Bool          -- ^ Result.
checkRepMapCompleteIMIS repMap imis =
  let -- All values, that appear in the union of all IntSet
      oldSets    = IS.unions imis
  in  checkRepMapCompleteIM repMap imis && checkRepMapCompleteIS repMap oldSets

{-|
Replace all 'IS.Key's from an 'IntSet' according to mappings from an 'IntMap' 'Int'. Entries, that
cannot be found in the 'IntMap' will no be changed.
-}
replaceIS ::
     IntMap Int -- ^ 'IntMap' containing the mapping from old 'IS.Key's to new 'IS.Key's.
  -> IntSet     -- ^ 'IntSet' to be modified.
  -> IntSet     -- ^ Resulting new 'IntSet'.
replaceIS repMap is =
  IS.map (\oK ->
    let nK = repMap !? oK
    in  fromMaybe oK nK
  ) is

{-|
Replace all 'IM.Key's in an 'IntMap' according to mappings from an 'IntMap' 'Int'. Entries that
cannot be found in the 'IntMap' will not be changed.
-}
replaceIMKeys ::
     IntMap Int -- ^ 'IntMap' containing the mapping from old 'IM.Key's to new 'IM.Key's.
  -> IntMap a   -- ^ 'IntMap' in which 'IM.Key's shall be replaced.
  -> IntMap a   -- ^ Resulting new 'IntMap' with replaced 'IM.Key's.
replaceIMKeys repMap im =
  IM.mapKeys (\oK ->
    let nK = repMap !? oK
    in  fromMaybe oK nK
  ) im

{-|
Replace all lookup-'IM.Key's of an 'IntMap' 'IntSet' and all values in the 'IntSet' by a given
mapping from an 'IntMap'. Entries that cannot be found in the 'IntMap' will not be changed.
-}
replaceIMIS ::
     IntMap Int    -- ^ 'IntMap' containing the mapping from old 'IS.Key's to new 'IS.Key's.
  -> IntMap IntSet -- ^ Original structure, which to replace both keys and values.
  -> IntMap IntSet -- ^ modified structure.
replaceIMIS repMap imis = undefined
  -- Replace all values in all IntSet
    IM.map (replaceIS repMap)
  -- Replace all lookup keys
  . replaceIMKeys repMap
  $ imis
