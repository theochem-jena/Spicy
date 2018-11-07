module Spicy.Trajectory
( filterByCriteria
, criterionDistance
, criterionAngle4Atoms
, findNearestAtom
) where
import           Control.Parallel.Strategies
import           Data.List
import           Data.Maybe
import           Lens.Micro.Platform
import           Spicy.Math
import           Spicy.Types
import           System.Random

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
