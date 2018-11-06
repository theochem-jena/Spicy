module Spicy.Trajectory
( filterByCriteria
, criterionDistance
) where
import           Control.Parallel.Strategies
import           Data.List
import           Data.Maybe
import           Spicy.Types
import           System.Random
import Spicy.Math
import           Lens.Micro.Platform

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
    joinBoolLists [] = []
    joinBoolLists [l] = l
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
criterionDistance :: (Int, Int) -> (Double -> Bool) -> Molecule -> Bool
criterionDistance (a, b) c m = c $ hmVecDistance (coord_a, coord_b)
  where
    coord_a = r3Vec2hmVec $ ((m ^. molecule_Atoms) !! a) ^. atom_Coordinates
    coord_b = r3Vec2hmVec $ ((m ^. molecule_Atoms) !! b) ^. atom_Coordinates
