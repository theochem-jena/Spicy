{-
This module deals with the partitioning of the system, creation of bonds,
assignment of substructeres to layers and creation of ghost atoms.
Following conventions shall apply:
  - highest level region has index 0
  - lowest level region has highest index
-}
module Spicy.MolecularSystem
( guessBonds
, insertPseudoBond
, isolateLayer
, filterByCriteria
, criterionDistance
, criterionAngle4Atoms
, findNearestAtom
) where
import           Control.Applicative
import           Data.List
import           Data.Map              (Map)
import qualified Data.Map              as Map
import           Data.Maybe
import           Lens.Micro.Platform
import qualified Numeric.LinearAlgebra as Algebra
import qualified Numeric.LinearAlgebra (C, Matrix, R, Vector)
import           Spicy.Math
import           Spicy.Types
import           System.Random
import           Control.Parallel.Strategies
import qualified Data.IntSet as I


--------------------------------------------------------------------------------
-- Manipulations of the Molecule
--------------------------------------------------------------------------------
-- | helper for the manipulation of bonds
data BondManipulation = Create | Delete deriving Eq

-- | Create or delete bonds between pairs of atoms and return a tuple with their
-- | updated connectivities.
-- | The atoms need to be given with their respective index in the molecule and
-- | the atom itseld
manipulateBondAtomic :: BondManipulation -> (Int, Atom) -> (Int, Atom) -> (Atom, Atom)
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

-- | Create or delete bonds between pairs of atoms and give back a molecule
manipulateBond :: BondManipulation -> Int -> Int -> Molecule -> Maybe Molecule
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

-- | For a molecule without connectivity yet, guess bonds based on interatomic
-- | distances only and give back the molecule with bonds
guessBonds :: Maybe Double -> Molecule -> Molecule
guessBonds scale mol = mol & molecule_Atoms .~ updatedAtoms
  where
    atoms = mol ^. molecule_Atoms
    atomIndRange = [ 0 .. length atoms - 1 ]
    updatedAtoms =
      [ (atoms !! a) & atom_Connectivity .~
      ( I.fromList . nub . concat $
        [ if isNothing $ cD (atoms !! a) (atoms !! b)
            then I.toList $ (atoms !! a) ^. atom_Connectivity
            else
              if  fromMaybe 1.4 scale * fromJust (cD (atoms !! a) (atoms !! b))
                  >= d (atoms !! a) (atoms !! b) && a/=b
                then b : (I.toList $ (atoms !! a) ^. atom_Connectivity)
                else I.toList $ (atoms !! a) ^. atom_Connectivity
        | b <- atomIndRange
        ]
      )
      | a <- atomIndRange
      ]
    d :: Atom -> Atom -> Double
    d a b =
      hmVecLength
        ( r3Vec2hmVec (a ^. atom_Coordinates)
        - r3Vec2hmVec (b ^. atom_Coordinates)
        )
    cD :: Atom -> Atom -> Maybe Double
    cD a b =
      (+) <$>
      Map.lookup (a ^. atom_Element) covalentRadii <*>
      Map.lookup (b ^. atom_Element) covalentRadii

-- | Insert a pseudo atom in a molecule with in a given scaled distance. Atom a
-- | and b are retained amd their bond will be untouched. A pseudo atom in
-- | between them will be inserted (at the end of the atom list of the molecule).
-- | The position of the pseudo atom (set 2 according to the paper) is
-- | calculated as described in https://doi.org/10.1016/S0166-1280(98)00475-8
insertPseudoBond :: Maybe Element -> Maybe Double -> Int -> Int -> Molecule -> Maybe Molecule
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

-- | Isolate parts of a molecule (defined by its indices) as a new ONIOM layer
-- | and saturate resulting dangling (single) bonds with a optionally specified
-- | element in a optionally specified distance.
-- | This if the fundamental principle for mechanical embedding. The pseudo
-- | get marked as such and can therefore be fixed in optimizations
-- | https://www.sciencedirect.com/science/article/pii/S0166128098004758?via%3Dihub
isolateLayer :: [Int] -> Maybe Element -> Maybe Double -> Molecule -> Maybe Molecule
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

{-
-- | Splits a molecule with correct bonds in individual fragments
fragmentMolecule :: Molecule -> SuperMolecule
fragmentMolecule m =
  where
    atoms = m ^. molecule_Atoms
    atomsBonds = map (^. atom_Connectivity) atoms
-}


--------------------------------------------------------------------------------
-- Analysis and filtering of molecules
--------------------------------------------------------------------------------
-- | Filter a trajectory by multiple criteria
filterByCriteria :: [Molecule -> Bool] -> Trajectory -> Trajectory
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

-- | Distance criterion (larger, smaller, equal, ...) of two atoms in a molecule
-- | to be fullfilled
criterionDistance :: (Int, Int) -> (Double -> Bool) -> Molecule -> Maybe Bool
criterionDistance (a, b) c m
  | a < nAtoms && b < nAtoms = Just $ c $ hmVecDistance (aCoord, bCoord)
  | otherwise = Nothing
  where
    atoms = m ^. molecule_Atoms
    nAtoms = length atoms
    aCoord = r3Vec2hmVec $ (atoms !! a) ^. atom_Coordinates
    bCoord = r3Vec2hmVec $ (atoms !! b) ^. atom_Coordinates

-- | Angle criterion, defined between 2x2 atoms. This is the general case, a2
-- | and b1 can be the same to have a 3 atom angle
criterionAngle4Atoms :: ((Int, Int), (Int, Int)) -> (Double -> Bool) -> Molecule -> Maybe Bool
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

-- | Give a point (might be coordinates of an atom) and find the closest atom to
-- | it. Gives a tuple of the distance to the atom, the index of the atom in the
-- | molecule and the atom itself
findNearestAtom :: R3Vec -> Molecule -> Maybe (Double, Int, Atom)
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

{-
fragmentMolecule :: Molecule -> SuperMolecule
fragmentMolecule m =
  where
    atoms = m ^. molecule_Atoms
    atomsBonds = map (^. atom_Connectivity) atoms
    bondPartners =
-}
--------------------------------------------------------------------------------
-- Generic Helper Functions
--------------------------------------------------------------------------------
-- | given a list of substituions [(new, old)] replace all elements in a list
-- | with the new ones and remove them if there is no correponding new element
substituteElemsInList :: Eq a => [(a, a)] -> [a] -> [a]
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

-- | Give a list of interesting indices for a list which is meant to be a subset
-- | of the original list and return a list where old index, new index and the
-- | retained elements are stored
remapIndices :: [Int] -> [a] -> Maybe [(Int, Int, a)]
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

-- | replace the Nth element from a list with a new element
replaceNth :: Int -> a -> [a] -> [a]
replaceNth n newElement oldList =
  take n oldList ++ [newElement] ++ drop (n + 1) oldList

-- | Check if two lists have some overlap
hasOverlap :: (Eq a) => [a] -> [a] -> Bool
hasOverlap a b = length (intersect a b) > 0


--------------------------------------------------------------------------------
-- tabulated data
--------------------------------------------------------------------------------
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
