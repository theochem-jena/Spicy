{-|
Module      : Spicy.MolecularSystem
Description : Handling molecular informations
Copyright   : Phillip Seeber, 2019
License     : GPL-3
Maintainer  : phillip.seeber@uni-jena.de
Stability   : experimental
Portability : POSIX, Windows

This module deals with the partitioning of the system, creation of bonds,
assignment of substructeres to layers and creation of ghost atoms.
Following conventions shall apply:

    * highest level region has highest index

    * lowest level region has index 0 and contains the complete system
-}
{-# LANGUAGE BangPatterns #-}
module Spicy.MolecularSystem
( -- * Manipulation of the Molecule
  guessBonds
, isolateLayer
, wrapFragmentsToBox
, ReplicationAxis(..)
, replicateSystemAlongAxis
, fragmentMolecule
, shiftFragment
  -- * Filtering trajectories and analysing molecules
, filterByCriteria
, criterionDistance
, criterionAngle4Atoms
, findNearestAtom
, FragmentBonds(..)
-- * Utility and Analysis functions
, distanceMatrix
) where
import           Control.Parallel.Strategies
import qualified Data.Array.IArray           as A
import           Data.IntSet                 (IntSet)
import qualified Data.IntSet                 as I
import           Data.List
import           Data.Map                    (Map)
import qualified Data.Map                    as Map
import           Data.Maybe
import           Lens.Micro.Platform
import qualified Numeric.LinearAlgebra       as Algebra
import           Spicy.Math
import           Spicy.Types


----------------------------------------------------------------------------------------------------
{-|
Auxiliary data type for the manipulation of bonds.
-}
data BondManipulation = Create | Delete deriving Eq

{-|
Create or delete bonds between pairs of 'Atom's A and B and return a tuple with their updated
connectivities.
-}
manipulateBondAtomic ::
     BondManipulation -- ^ Create or delete a bond.
  -> (Int, Atom)      -- ^ (Index of 'Atom' A in the 'Molecule', 'Atom' A).
  -> (Int, Atom)      -- ^ (Index of 'Atom' B in the 'Molecule', 'Atom' B).
  -> (Atom, Atom)     -- ^ (Updated 'Atom' A, updated 'Atom' B).
manipulateBondAtomic action (i_a, a) (i_b, b) = (atomNew_a, atomNew_b)
  where
    atom_a = a
    atom_b = b
    connectivity_a = I.toList $ atom_a ^. atom_Connectivity
    connectivity_b = I.toList $ atom_b ^. atom_Connectivity
    connectivityCleaned_a =
      if action == Create
        then i_b : connectivity_a
        else delete i_b connectivity_a
    connectivityCleaned_b =
      if action == Create
        then i_a : connectivity_b
        else delete i_a connectivity_b
    atomNew_a = atom_a & atom_Connectivity .~ I.fromList connectivityCleaned_a
    atomNew_b = atom_b & atom_Connectivity .~ I.fromList connectivityCleaned_b

{-|
Create or delete bonds between pairs of 'Atom's A and B and give back an updated 'Molecule'.
-}
manipulateBond ::
     BondManipulation -- ^ Create or delete a bond.
  -> Int              -- ^ Index of 'Atom' A.
  -> Int              -- ^ Index of 'Atom' B.
  -> Molecule         -- ^ The 'Molecule' which to update.
  -> Maybe Molecule   -- ^ Result 'Molecule'. Might fail with 'Nothing' if indices are invalid.
manipulateBond action i_a i_b mol
  | i_a > maxInd || i_a < 0 || i_b > maxInd || i_b < 0 = Nothing
  | otherwise = Just $ mol & molecule_Atoms .~ atomsNew
  where
    maxInd = length atoms - 1
    atoms = mol ^. molecule_Atoms
    atom_a = atoms !! i_a
    atom_b = atoms !! i_b
    (atomNew_a, atomNew_b) = manipulateBondAtomic action (i_a, atom_a) (i_b, atom_b)
    atomsNew = (replaceNth i_b atomNew_b) . (replaceNth i_a atomNew_a) $ atoms

{-|
For a 'Molecule' without connectivity yet, guess bonds based on interatomic distances only and give
back the 'Molecule' with bonds.
-}
guessBonds ::
     Maybe Double -- ^ An optional scaling factor for the covalent radii.
  -> Molecule     -- ^ 'Molecule', of which to guess the bonds.
  -> Molecule     -- ^ 'Molecule' with new bonds guessed.
guessBonds scale mol = mol & molecule_Atoms .~ updatedAtoms
  where
    !distances = distanceMatrix mol
    atoms = mol ^. molecule_Atoms
    atomIndRange = [ 0 .. length atoms - 1 ]
    indexedAtoms = zip atomIndRange atoms
    updatedAtoms = zipWith (\refAt otherAts -> makeAllAtomsBonds distances refAt otherAts)
      indexedAtoms (replicate (length atoms) indexedAtoms) `using` (parListChunk 500) rdeepseq
    -- This function takes a single reference atom (with index) and a indexed list of other atoms.
    --  It adds all bonds to other atoms found to the reference atom
    makeAllAtomsBonds :: A.Array (Int, Int) Double -> (Int, Atom) -> [(Int, Atom)] -> Atom
    makeAllAtomsBonds distMat thisAt allAt =
      snd $
      foldl (\origAtom otherAtom -> if (fst origAtom) == (fst otherAtom)
        then origAtom
        else if fromMaybe 1.4 scale * fromJust (cD (snd origAtom) (snd otherAtom)) <=
                (distMat A.! (fst origAtom, fst otherAtom))
               then origAtom
               else (fst origAtom, (snd origAtom) & atom_Connectivity %~ I.insert (fst otherAtom))
      ) thisAt allAt
    cD :: Atom -> Atom -> Maybe Double
    cD a b =
      (+) <$>
      Map.lookup (a ^. atom_Element) covalentRadii <*>
      Map.lookup (b ^. atom_Element) covalentRadii

{-|
Insert a pseudo 'Atom' in a 'Molecule' with in a given scaled distance. 'Atom's A and B are retained
and their bond will be untouched. A pseudo 'Atom' in between them will be inserted (at the end of
the 'Atom' list of the 'Molecule'). The position of the pseudo atom (set 2 according to the paper)
is calculated as described in <https://doi.org/10.1016/S0166-1280(98)00475-8>.
-}
insertPseudoBond ::
     Maybe Element   -- ^ Optionally a chemical 'Element' to be used as pseubond capping. Hydrogen
                     --   is used if not specified.
  -> Maybe Double    -- ^ Optional scaling factor, where to place the pseudoatom on the original
                     --   bond vector. 1 means the pseudoatom will be placed, where the original cut
                     --   off atom was, 0.5 by half the distance of the original bond. Default is to
                     --   calculate the scaling by covalent radii of the atoms.
  -> Int             -- ^ Index of 'Atom' A.
  -> Int             -- ^ Index of 'Atom' B.
  -> Molecule        -- ^ Original 'Molecule' to modify.
  -> Maybe Molecule  -- ^ Result 'Molecule' with pseudo 'Atom' and pseubond inserted. May fail
                     --   because of out of bound indices or because covalent radii for the elements
                     --   are unknown.
insertPseudoBond psElementTemplate psScalingTemplate i_a i_b mol
  | i_a > maxInd || i_a < 0 || i_b > maxInd || i_b < 0 = Nothing
  | (isNothing covRad_ab || isNothing covRad_ap) && isNothing psElementTemplate = Nothing
  | isNothing (hmVec2r3Vec r_p) = Nothing
  | otherwise = Just $ mol & molecule_Atoms .~ (atoms ++ [atom_p])
  where
    maxInd = length atoms - 1
    atoms = mol ^. molecule_Atoms
    atom_a = atoms !! i_a
    atom_b = atoms !! i_b
    -- select an element to saturate the dangling bond, defaulting to hydrogen
    psElement = fromMaybe H psElementTemplate
    -- pseudo atom will be placed at r_p = r_a + g * (r_b - r_a)
    -- g (psScaling) is either given or calculated from the covalent radii of
    -- a, b and p
    -- atom a is the inner layer, p the pseudo atom and b the outer layer
    covRad_p = Map.lookup psElement covalentRadii
    covRad_a = Map.lookup (atom_a ^. atom_Element) covalentRadii
    covRad_b = Map.lookup (atom_b ^. atom_Element) covalentRadii
    covRad_ab = (+) <$> covRad_a <*> covRad_b
    covRad_ap = (+) <$> covRad_a <*> covRad_p
    psScaling = fromMaybe (fromJust covRad_ap / fromJust covRad_ab) psScalingTemplate -- this should be safe as i check covRad in the guards
    r_a = r3Vec2hmVec $ atom_a ^. atom_Coordinates
    r_b = r3Vec2hmVec $ atom_b ^. atom_Coordinates
    r_ab = r_b - r_a
    r_p = r_a + Algebra.vector [psScaling] * r_ab
    atom_p = Atom
      { _atom_Element = psElement
      , _atom_Label = ""
      , _atom_IsPseudo = True
      , _atom_FFType = ""
      , _atom_PCharge = Nothing
      , _atom_Coordinates = fromJust $ hmVec2r3Vec r_p
      , _atom_Connectivity = I.fromList [i_a]
      }

{-|
Isolate parts of a molecule (defined by its indices) as a new ONIOM layer and saturate resulting
dangling (single) bonds with a optionally specified element in a optionally specified distance. This
is the fundamental principle for mechanical embedding. The pseudo get marked as such and can
therefore be fixed in optimisations
<https://www.sciencedirect.com/science/article/pii/S0166128098004758?via%3Dihub>.
-}
isolateLayer ::
     [Int]          -- ^ Indices of the 'Atom's in the new ONIOM layer.
  -> Maybe Element  -- ^ Optional element to be used as capping atom. See 'insertPseudoBond'
  -> Maybe Double   -- ^ Optional scaling factor for pseudobonds. See 'insertPseudoBond'.
  -> Molecule       -- ^ Molecule from which to cut the ONIOM layer.
  -> Maybe Molecule -- ^ New pseudoatom saturated ONIOM layer molecule.
isolateLayer nlInd pseudoElement pseudoScaling mol
  | maximum nlInd > maximum olIndRange || minimum nlInd < 0 = Nothing
  | Nothing `elem` pseudoAtomsToOl = Nothing
  | isNothing replacementList = Nothing
  | otherwise = Just nlMolecule
  where
    olAtoms = mol ^. molecule_Atoms
    olIndRange = [0 .. length olAtoms - 1 ]
    nlAtoms = [ olAtoms !! i | i <- nlInd ] :: [Atom]
    olAtomicBonds = map (I.toList . (^. atom_Connectivity)) olAtoms             -- atomwise list of bonds to other atoms of the old layer
    nlAtomicBonds = map (I.toList . (^. atom_Connectivity)) nlAtoms
    pseudoAtomsToOl =                                                           -- the set of pseudo atoms that are added to the outer layer when separating the inner layers
      concat
      [ [ last . _molecule_Atoms <$>                                            -- take the last atom added to the molecule
          insertPseudoBond pseudoElement pseudoScaling nA roA mol               -- which is added for a pseudo bond to all new layer atoms
        | roA <-
          foldl (\rb b -> if b `elem` nlInd
                            then rb
                            else b:rb)
          [] (olAtomicBonds !! nA)                                              -- if it connects to any old layer atom, which is not a new layer atom
        ]
      | nA <- nlInd
      ] :: [Maybe Atom]
    -- old layer with pseudo atoms added
    pseudoMolAtoms = olAtoms ++ map fromJust pseudoAtomsToOl
    -- now start removing the old layer and get informations how the index remapping is
    replacementList =
      remapIndices
      (nlInd ++ [ last olIndRange + 1 .. length pseudoMolAtoms - 1 ])
      pseudoMolAtoms
    -- only the new set of atoms but without the updated connectivities
    newAtomsOldConnectivity = map (\(a, b, c) -> c) $ fromJust replacementList
    newAtomsOldConnectivitiesOnly = map (I.toList . (^. atom_Connectivity)) newAtomsOldConnectivity
    -- the replacment list for updating connectivities
    indMappingList = map (\(a, b, c) -> (a, b)) $ fromJust replacementList
    newAtomsNewConnectivitiesOnly =
      map (substituteElemsInList indMappingList) newAtomsOldConnectivitiesOnly
    -- the new atoms with updated connectivities
    newAtoms =
      [ (newAtomsOldConnectivity !! i) & atom_Connectivity .~ I.fromList (newAtomsNewConnectivitiesOnly !! i)
      | i <- [0 .. length newAtomsOldConnectivity - 1 ]
      ]
    nlMolecule = mol & molecule_Atoms .~ newAtoms

{-|
Wrap fragments moleculewise to the rectangular unit cell. This function assumes that the origin is
at (0, 0, 0).
-}
wrapFragmentsToBox ::
     R3Vec         -- ^ Length of x y and z axes of a rectangular unit cell.
  -> SuperMolecule -- ^ Fragmented 'Supermolecule' before wrapping.
  -> SuperMolecule -- ^ Wrapped 'Supermolecule' after wrapping.
wrapFragmentsToBox (x, y, z) (supermol, fragments) = (updatedSupermol, wrappedFragments)
  where
    wrappedFragments = map (shiftFragmentToUnitCell (x, y, z)) fragments
    fragmentAtoms = map (^. molecule_Atoms) wrappedFragments
    updatedSupermol = supermol & molecule_Atoms .~ concat fragmentAtoms

{-|
Select the axis along which to replicate the system.
-}
data ReplicationAxis = AxisX | AxisY | AxisZ deriving Eq

{-|
Replicate along a given Axis. Shift coordinates, so that the new unit cell has its origin allways at
0,0,0.
-}
replicateSystemAlongAxis ::
     R3Vec           -- ^ Length of x y and z axes of a rectangular unit cell.
  -> ReplicationAxis -- ^ Axis along which to replicate the system.
  -> Molecule        -- ^ 'Molecule' to replicate.
  -> Molecule        -- ^ Replication result 'Molecule'.
replicateSystemAlongAxis (bx, by, bz) axis m = mReplicatedShifted
  where
    atoms = m ^. molecule_Atoms
    nAtoms = length atoms
    replicaAtomsPos =
      [ a
        & atom_Coordinates .~ case axis of
          AxisX -> fromJust . hmVec2r3Vec $ (r3Vec2hmVec $ a ^. atom_Coordinates) + (r3Vec2hmVec (bx, 0, 0))
          AxisY -> fromJust . hmVec2r3Vec $ (r3Vec2hmVec $ a ^. atom_Coordinates) + (r3Vec2hmVec (0, by, 0))
          AxisZ -> fromJust . hmVec2r3Vec $ (r3Vec2hmVec $ a ^. atom_Coordinates) + (r3Vec2hmVec (0, 0, bz))
        & atom_Connectivity .~ (I.fromList . (map (+ nAtoms)) $ (I.toList $ a ^. atom_Connectivity))
      | a <- atoms
      ]
    replicaAtomsNeg =
      [ a
        & atom_Coordinates .~ case axis of
          AxisX -> fromJust . hmVec2r3Vec $ (r3Vec2hmVec $ a ^. atom_Coordinates) - (r3Vec2hmVec (bx, 0, 0))
          AxisY -> fromJust . hmVec2r3Vec $ (r3Vec2hmVec $ a ^. atom_Coordinates) - (r3Vec2hmVec (0, by, 0))
          AxisZ -> fromJust . hmVec2r3Vec $ (r3Vec2hmVec $ a ^. atom_Coordinates) - (r3Vec2hmVec (0, 0, bz))
        & atom_Connectivity .~ (I.fromList . (map (+ (2 * nAtoms))) $ (I.toList $ a ^. atom_Connectivity))
      | a <- atoms
      ]
    allNewAtoms = atoms ++ replicaAtomsPos ++ replicaAtomsNeg
    reShiftVec =
      r3Vec2hmVec $
      case axis of
        AxisX -> (bx, 0, 0)
        AxisY -> (0, by, 0)
        AxisZ -> (0, 0, bz)
    allNewAtomsNewPos =
      [ a & atom_Coordinates .~ fromJust (hmVec2r3Vec (r3Vec2hmVec (a ^. atom_Coordinates) + reShiftVec))
      | a <- allNewAtoms
      ]
    mReplicatedShifted = m & molecule_Atoms .~ allNewAtomsNewPos


----------------------------------------------------------------------------------------------------
-- Analysis and filtering of molecules

{-|
Filter a 'Trajectory' by multiple criteria. All crieteria must be fullfilled for a 'Molecule'/frame
to be in the output.
-}
filterByCriteria ::
     [Molecule -> Bool] -- ^ A list of functions, that decide if a 'Molecule' will be kept or
                        --   discarded a.k.a. the criteria.
  -> Trajectory         -- ^ Input 'Trajectory' to filter.
  -> Trajectory         -- ^ Filtered output 'Trajectory'.
filterByCriteria cs t = resultTraj
  where
    -- Take criterion (c) and a list (l) and check for every element in the list
    -- in parallel if the criterion is fullfilled.
    filter2Bool :: (a -> Bool) -> [a] -> [Bool]
    filter2Bool c l = map c l `using` parList rdeepseq

    -- Filter lists are joined by logical AND
    joinBoolLists :: [[Bool]] -> [Bool]
    joinBoolLists []     = []
    joinBoolLists [l]    = l
    joinBoolLists (l:ls) = zipWith (&&) l (joinBoolLists ls)

    -- Boolean lists of the frames by criteria
    individualBoolLists = map (\e -> filter2Bool e t) cs `using` parList rdeepseq
    boolList = joinBoolLists individualBoolLists

    -- Mark each frame for deletion or keeping
    markedTraj = zip boolList t

    -- Filtered marked trajectory
    resultTraj = map snd . filter (\(b, m) -> b == True) $ markedTraj

{-|
Distance criterion (larger, smaller, equal, ...) of two atoms in a molecule.
-}
criterionDistance ::
     (Int, Int)       -- ^ Indices of the 2 'Atom's in the molecule, whos distance is calculated.
  -> (Double -> Bool) -- ^ A distance criterion, e.g. @ (> 5.0) @.
  -> Molecule         -- ^ 'Molecule' to check.
  -> Maybe Bool       -- ^ Filter result.
criterionDistance (a, b) c m
  | a < nAtoms && b < nAtoms = Just $ c $ hmVecDistance (aCoord, bCoord)
  | otherwise = Nothing
  where
    atoms = m ^. molecule_Atoms
    nAtoms = length atoms
    aCoord = r3Vec2hmVec $ (atoms !! a) ^. atom_Coordinates
    bCoord = r3Vec2hmVec $ (atoms !! b) ^. atom_Coordinates

{-|
Angle criterion, defined between 2x2 a1, a2 and b1 and b2 atoms. This is the general case, a2 and b1
can be the same to have a 3 atom angle.
-}
criterionAngle4Atoms ::
     ((Int, Int), (Int, Int)) -- ^ Index of the 'Atom's ((a1, a2), (b1, b2)).
  -> (Double -> Bool)         -- ^ An angle criterion, e.g. @ (< 0.5 * pi) @.
  -> Molecule                 -- ^ 'Molecule' to check.
  -> Maybe Bool               -- ^ Filter result.
criterionAngle4Atoms ((a1, a2), (b1, b2)) c m
  | a1 < nAtoms && a2 < nAtoms && b1 < nAtoms && b2 < nAtoms = Just $ c $ hmVecAngle (aVec, bVec)
  | otherwise = Nothing
  where
    atoms = m ^. molecule_Atoms
    nAtoms = length atoms
    a1Coord = r3Vec2hmVec $ (atoms !! a1) ^. atom_Coordinates
    a2Coord = r3Vec2hmVec $ (atoms !! a2) ^. atom_Coordinates
    b1Coord = r3Vec2hmVec $ (atoms !! b1) ^. atom_Coordinates
    b2Coord = r3Vec2hmVec $ (atoms !! b2) ^. atom_Coordinates
    aVec = a2Coord - a1Coord
    bVec = b2Coord - b1Coord

{-|
Give a point (might be coordinates of an 'Atom') and find the closest 'Atom' to it. Gives a tuple of
the distance to the 'Atom', the index of the 'Atom' in the 'Molecule' and the 'Atom' itself.
-}
findNearestAtom ::
     R3Vec                     -- ^ Position in space from which to search the nearest 'Atom'.
  -> Molecule                  -- ^ Molecule in which to look for the nearest 'Atom'.
  -> Maybe (Double, Int, Atom) -- ^ Result with (distance to the nearest 'Atom', index of the
                               --   nearest 'Atom', nearest 'Atom')
findNearestAtom pos m
  | length (m ^. molecule_Atoms) < 1 = Nothing
  | isNothing indexMaybe = Nothing
  | otherwise = Just (distances !! index, index, atoms !! index)
  where
    atoms = m ^. molecule_Atoms
    nAtoms = length atoms
    posVec = r3Vec2hmVec pos
    distances =
      [ hmVecDistance (posVec, r3Vec2hmVec $ i ^. atom_Coordinates)
      | i <- atoms
      ]
    smallestDistance = minimum distances
    indexMaybe = findIndex (== smallestDistance) distances
    index = fromJust indexMaybe

{-|
When segmenting a 'Molecule' into fragments, how should the bonds in the fragment be handled?
 SuperAndFragment -> In supermolecule the original bonds remain but in the fragments all bonds
   are remapped to have intrafragment bonds as they were in the super molecule
 NewGuess -> Applies bond guessing based on covalent radii fragmentwise only in the fragments
 RemoveAll -> Remove all bonds from supermolecule and fragments
 KeepBonds -> Do not change the bonds in the fragment and keep the indices from the supermol
-}
data FragmentBonds =
    OnlySuper               -- ^ Remove all bonds in the fragments and only keep them in the
                            --   'SuperMolecule'
  | SuperAndFragment        -- ^ In 'Supermolecule' the original bonds remain but in the fragments
                            --   all bonds are remapped to have intrafragment bonds as they were in
                            --   the 'SuperMolecule'.
  | NewGuess (Maybe Double) -- ^ Applies bond guessing based on covalent radii fragmentwise only in
                            --   the fragments.
  | RemoveAll               -- ^ Remove all bonds from 'Supermolecule' and fragments.
  | KeepBonds               -- ^ Do not change the bonds in the fragment and keep the indices from
                            --   the 'SuperMolecule'.
  deriving Eq

{-|
Detects fragment based on bond analysis. Bonds are handled according to 'FragmentBonds'.
-}
fragmentMolecule ::
     FragmentBonds       -- ^ Mode to handle bonds when fragmenting.
  -> Molecule            -- ^ 'Molecule' to fragment.
  -> Maybe SuperMolecule -- ^ Fragmented result 'SuperMolecule'.
fragmentMolecule bondHandling m
  | bondHandling == SuperAndFragment = Nothing
  | isNothing fragmentsBondsUpdate = Nothing
  | bondHandling == RemoveAll = Just (moleculeWithoutBonds, fromJust fragmentsBondsUpdate)
  | otherwise = Just (m, fromJust fragmentsBondsUpdate)
  where
    atoms = m ^. molecule_Atoms
    nAtoms = length atoms
    atomsBonds = map (^. atom_Connectivity) atoms
    bondPartners =
      [ I.insert i (atomsBonds !! i)
      | i <- [0 .. nAtoms - 1]
      ]
    fragmentsIndices = reduceToZeroOverlap bondPartners
    fragments =
      [ Molecule
          { _molecule_Label    = "Fragment " ++ show i
          , _molecule_Atoms    = [(m ^. molecule_Atoms) !! a | a <- I.toList (fragmentsIndices !! i)]
          , _molecule_Energy   = Nothing
          , _molecule_Gradient = Nothing
          , _molecule_Hessian  = Nothing
          }
      | i <- [0 .. length fragmentsIndices - 1]
      ]
    fragmentsBondsUpdate = case bondHandling of
      OnlySuper ->
        Just fragmentsWithoutBonds
      SuperAndFragment ->
        Nothing
      NewGuess s ->
        Just
        [ guessBonds s f
        | f <- fragments
        ]
      RemoveAll ->
        Just fragmentsWithoutBonds
      KeepBonds ->
        Just fragments
      where
        fragmentsWithoutBonds =
          [ f & molecule_Atoms .~
            map (\a -> a & atom_Connectivity .~ I.empty) (f ^. molecule_Atoms)
          | f <- fragments
          ]
    moleculeWithoutBonds =
      m & molecule_Atoms .~
      map (\a -> a & atom_Connectivity .~ I.empty) (m ^. molecule_Atoms)


{-|
Calculate the distance matrix of a 'Molecule'.
-}
distanceMatrix ::
     Molecule                  -- ^ Input 'Molecule'.
  -> A.Array (Int, Int) Double -- ^ Distance matrix of the 'Molecule'.
distanceMatrix mol = A.array
  ((0, 0), (length coordinates - 1, length coordinates - 1)) $
  rightUpperEntries ++ diagonalEntries ++ leftLowerEntries
  where
    coordinates = map (r3Vec2hmVec . (^. atom_Coordinates)) $ mol ^. molecule_Atoms
    rightUpperIndices =
      [ (i, j)
      | i <- [0 .. length coordinates - 1], j <- [0 .. length coordinates - 1], j > i
      ]
    !rightUpperEntries =
      map (\(i, j) ->
        ( (i, j)
        , hmVecDistance (coordinates !! i, coordinates !! j)
        )
      ) rightUpperIndices -- `using` (parListChunk 5000) rdeepseq
    diagonalEntries = [ ((i, i), 0.0) | i <- [0 .. length coordinates - 1]]
    leftLowerIndices = map (\(a, b) -> (b, a)) rightUpperIndices -- `using` (parListChunk chunkSize) rdeepseq
    !leftLowerEntries =
      zipWith (\(_, dist) inds -> (inds, dist)) rightUpperEntries leftLowerIndices
      -- `using` (parListChunk 5000) rdeepseq


----------------------------------------------------------------------------------------------------
-- Generic Helper Functions

{-|
For a given base 'Vector' defining a rectangular cell with the coordinate origin and a fragment
('Molecule') shift the 'Molecule' such, that it is at least partially contained in the unit cell.
-}
shiftFragmentToUnitCell ::
     R3Vec    -- ^ Cell 'Vector' of a rectangular unit cell.
  -> Molecule -- ^ Fragment that shall be shifted to the unit cell.
  -> Molecule -- ^ Shifted fragment.
shiftFragmentToUnitCell (bx, by, bz) m = moleculeShifted
  where
    atoms = m ^. molecule_Atoms
    atomCoords = map (^. atom_Coordinates) atoms
    xShifts =
      [ floor ((a ^. _1) / bx)
      | a <- atomCoords
      ] :: [Int]
    yShifts =
      [ floor ((a ^. _2) / by)
      | a <- atomCoords
      ] :: [Int]
    zShifts =
      [ floor ((a ^. _3) / bz)
      | a <- atomCoords
      ] :: [Int]
    extremaX = (minimum xShifts, maximum xShifts)
    extremaY = (minimum yShifts, maximum yShifts)
    extremaZ = (minimum zShifts, maximum zShifts)
    effectiveX =
      if (0 `elem` xShifts)
        then 0
        else if (abs . fst $ extremaX) > (snd extremaX)
               then fst extremaX
               else snd extremaX
    effectiveY =
      if (0 `elem` yShifts)
        then 0
        else if (abs . fst $ extremaY) > (snd extremaY)
               then fst extremaY
               else snd extremaY
    effectiveZ =
      if (0 `elem` zShifts)
        then 0
        else if (abs . fst $ extremaZ) > (snd extremaZ)
               then fst extremaZ
               else snd extremaZ
    xShift = (-1) * (fromIntegral effectiveX) * bx
    yShift = (-1) * (fromIntegral effectiveY) * by
    zShift = (-1) * (fromIntegral effectiveZ) * bz
    shiftVec = (xShift, yShift, zShift)
    moleculeShifted = shiftFragment shiftVec m

{-|
Check for a fragment if if is partially in the unit cell (one atom in the unit cell is enough).
-}
isFragInUnitCell :: R3Vec -> Molecule -> Bool
isFragInUnitCell (bx, by, bz) m = True `elem` atomsInUnitCell
  where
    atoms = m ^. molecule_Atoms
    atomCoords = map (^. atom_Coordinates) atoms
    atomsInUnitCell =
      [ (\(ax, ay, az) ->
          (ax / bx) >= 0 && (ax / bx) <= 1.0 &&
          (ay / by) >= 0 && (ay / by) <= 1.0 &&
          (az / bz) >= 0 && (az / bz) <= 1.0
        ) (a ^. atom_Coordinates)
      | a <- atoms
      ]

{-|
Shift an 'R3Vec' back to the boundaries of a rectangular box.
-}
toBaseVec ::
    R3Vec -- ^ Point in space, possibly outside a unit cell.
 -> R3Vec -- ^ Base vector of the cell.
 -> R3Vec -- ^ Point shifted back to the first unit cell.
toBaseVec (cx, cy, cz) (bx, by, bz) = coordR3Shifted
  where
    -- in each axis, how many unit vectors offset could i shift the molecule back
    (ox, oy, oz) =
      ( floor (cx / bx)
      , floor (cy / by)
      , floor (cz / bz)
      )
    coordR3Shifted =
      ( cx + (fromIntegral ox) * bx
      , cy + (fromIntegral oy) * by
      , cz + (fromIntegral oz) * bz
      )

{-|
Shift a 'Molecule' by a given vector. The shift only applies to the cartesian coordinates.
-}
shiftFragment ::
     R3Vec    -- ^ Vector to shift.
  -> Molecule -- ^ 'Molecule' to shift.
  -> Molecule -- ^ Shifted 'Molecule'.
shiftFragment shiftVec m = moleculeShifted
  where
    atoms = m ^. molecule_Atoms
    atomCoords = map (^. atom_Coordinates) atoms
    atomCoordsHM = map r3Vec2hmVec atomCoords
    shiftHMVec = r3Vec2hmVec shiftVec
    atomCoordsShiftedHM = map (+ shiftHMVec) atomCoordsHM
    atomCoordsShifted = map (fromJust . hmVec2r3Vec) atomCoordsShiftedHM
    atomsShifted =
      [ (atoms !! i) & atom_Coordinates .~ (atomCoordsShifted !! i)
      | i <- [0 .. length atoms - 1]
      ]
    moleculeShifted = m & molecule_Atoms .~ atomsShifted

{-|
Given a list of substituions @[(new, old)]@ replace all elements in a list with the new ones and
remove them if there is no correponding new element.
-}
substituteElemsInList ::
     Eq a
  => [(a, a)] -- ^ Substitution pairs in the form @(new, old)@.
  -> [a]      -- ^ Original list to substitute.
  -> [a]      -- ^ Result list with substituions applied.
substituteElemsInList substiList origList = substitutedList
  where
    new = map fst substiList
    old = map snd substiList
    intermediateSubstitutedList =
      foldr (\x acc ->
        if x `elem` old
          then
            (new !! fromJust (elemIndex x old), old !! fromJust (elemIndex x old)) : acc
          else acc
        )
        [] origList
    substitutedList = map fst intermediateSubstitutedList

{-|
Give a list of interesting indices for a list which is meant to be a subset of the original list and
return a list where old index, new index and the retained elements are stored.
-}
remapIndices ::
     [Int]                 -- ^ List of indices for elements from the old list to keep.
  -> [a]                   -- ^ The old list.
  -> Maybe [(Int, Int, a)] -- ^ New List with @(newIndex, oldIndex, element)@.
remapIndices ind origList
  | maximum ind > (length origList - 1) || minimum ind < 0 = Nothing
  | otherwise = Just newIndOrigIndNewElemList
  where
    origIndNewElemList = [ (i, origList !! i) | i <- ind ]
    newIndList = [ 0 .. length origIndNewElemList - 1 ]
    newIndOrigIndNewElemList =
      [ (nI, fst (origIndNewElemList !! nI), snd (origIndNewElemList !! nI))
      | nI <- newIndList
      ]

{-|
Replace the Nth element from a list with a new element.
-}
replaceNth ::
     Int -- ^ The Nth element to replace.
  -> a   -- ^ The substituion element.
  -> [a] -- ^ Original list.
  -> [a] -- ^ Replaced list.
replaceNth n newElement oldList =
  take n oldList ++ [newElement] ++ drop (n + 1) oldList

{-|
Delete the Nth element of a list.
-}
deleteNth :: Int -> [a] -> [a]
deleteNth n l = (take n l) ++ (drop (n + 1) l)

{-|
This reduces a list of IntSets, that potentially has some overlap in the set elements, to a list of
sets, that has no elements in common, by recursively combining all sets, that have overlap.
-}
reduceToZeroOverlap :: [IntSet] -> [IntSet]
reduceToZeroOverlap [] = []
reduceToZeroOverlap [a] = [a]
reduceToZeroOverlap (a:as) =
  if (all null overlaps)
    then (a:as)
    else reduceToZeroOverlap reducedSet
  where
    -- From a list of possible sets, which with the first could have overlap, the results will be
    -- all, that have an overlap
    setsWithOverlap :: IntSet -> [IntSet] -> [IntSet]
    setsWithOverlap b bl = [ i | i <- bl, not . I.null $ I.intersection b i, b /= i]

    -- List of lists of sets with overlap to a fragment
    -- (fragments -> other fragments that have overlap -> atoms in this fragment)
    overlaps =
      [ setsWithOverlap ((a:as) !! i) (deleteNth i (a:as))
      | i <- [0 .. length (a:as) - 1]
      ] :: [[IntSet]]

    -- Join all sets with overlap at once
    reducedSet =
      nub
      [ I.unions (((a:as) !! i) : (overlaps !! i))
      | i <- [0 .. length (a:as) - 1]
      ]


----------------------------------------------------------------------------------------------------
-- tabulated data
{-|
The covalent radii of alle elements in Angstrom.
-}
covalentRadii :: Map Element Double
covalentRadii = Map.fromList
  [ (H  , 0.31 ),                                                                                                                                                                                                                                 (He , 0.28 ),
    (Li , 1.28 ), (Be , 0.96 ),                                                                                                                                             (B  , 0.84 ), (C  , 0.76 ), (N  , 0.71 ), (O  , 0.66 ), (F  , 0.57 ), (Ne , 0.58 ),
    (Na , 1.66 ), (Mg , 1.41 ),                                                                                                                                             (Al , 1.21 ), (Si , 1.11 ), (P  , 1.07 ), (S  , 1.05 ), (Cl , 1.02 ), (Ar , 1.06 ),
    (K  , 2.03 ), (Ca , 1.76 ), (Sc , 1.70 ), (Ti , 1.60 ), (V  , 1.53 ), (Cr , 1.39 ), (Mn , 1.61 ), (Fe , 1.52 ), (Co , 1.50 ), (Ni , 1.24 ), (Cu , 1.32 ), (Zn , 1.22 ), (Ga , 1.22 ), (Ge , 1.20 ), (As , 1.19 ), (Se , 1.20 ), (Br , 1.20 ), (Kr , 1.16 ),
    (Rb , 2.20 ), (Sr , 1.95 ), (Y  , 1.90 ), (Zr , 1.75 ), (Nb , 1.64 ), (Mo , 1.54 ), (Tc , 1.47 ), (Ru , 1.46 ), (Rh , 1.42 ), (Pd , 1.39 ), (Ag , 1.45 ), (Cd , 1.44 ), (In , 1.42 ), (Sn , 1.39 ), (Sb , 1.39 ), (Te , 1.38 ), (I  , 1.39 ), (Xe , 1.40 ),
    (Cs , 2.44 ), (Ba , 2.15 ), (Lu , 1.87 ), (Hf , 1.75 ), (Ta , 1.70 ), (W  , 1.62 ), (Re , 1.51 ), (Os , 1.44 ), (Ir , 1.41 ), (Pt , 1.36 ), (Au , 1.36 ), (Hg , 1.32 ), (Tl , 1.45 ), (Pb , 1.46 ), (Bi , 1.48 ), (Po , 1.40 ), (At , 1.50 ), (Rn , 1.50 ),
    (Fr , 2.60 ), (Ra , 2.21 ),

    (La , 2.07 ), (Ce , 2.04 ), (Pr , 2.03 ), (Nd , 2.01 ), (Pm , 1.99 ), (Sm , 1.98 ), (Eu , 1.98 ), (Gd , 1.96 ), (Tb , 1.50 ), (Dy , 1.92 ), (Ho , 1.92 ), (Er , 1.89 ), (Tm , 1.90 ), (Yb , 1.87 ),
    (Ac , 1.92 ), (Th , 2.06 ), (Pa , 2.00 ), (U  , 1.96 ), (Np , 1.90 ), (Pu , 1.87 ), (Am , 1.80 ), (Cm , 1.69 )
  ]
