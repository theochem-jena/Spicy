{-|
Module      : Spicy.MolWriter
Description : Writing chemical data formats
Copyright   : Phillip Seeber, 2019
License     : GPL-3
Maintainer  : phillip.seeber@uni-jena.de
Stability   : experimental
Portability : POSIX, Windows

A module which converts the internal Molecule representation to a string, which is a common chemical
file format, that can be read by Avogadro, VMD, OpenBabel etc.. The writers are not fool proof with
respect to force field types, which should always be remembered when usings its results.
-}
{-# LANGUAGE ScopedTypeVariables #-}
module Spicy.Writer.Molecule
( writeXYZ
, writeTXYZ
, writeMOL2
, writePDB
, writeSpicy
) where
import           Data.Aeson.Encode.Pretty
import           Data.Attoparsec.Text.Lazy (isEndOfLine)
import           Data.Either
import           Data.Foldable
import qualified Data.IntMap.Lazy          as IM
import qualified Data.IntSet               as IS
import           Data.Maybe
import           Data.Sequence             (Seq)
import qualified Data.Sequence             as S
import           Data.Text.Lazy            (Text)
import qualified Data.Text.Lazy            as T
import qualified Data.Text.Lazy.Encoding   as T
import           Lens.Micro.Platform
import           Prelude                   hiding (cycle, foldl1, foldr1, head,
                                            init, last, maximum, minimum, tail,
                                            take, takeWhile, (!!))
import           Spicy.Molecule.Util
import           Spicy.Types
import           Text.Printf


{-|
Write a Molden XYZ file from a molecule. This format ignores all deep level layers of a molecule.
-}
writeXYZ :: Molecule -> Either String Text
writeXYZ mol = toXYZ <$> checkMolecule mol
  where
    -- Assemble a XYZ file from a molecule
    toXYZ :: Molecule -> Text
    toXYZ m =
      let -- Line 1 of the header contains the number of atoms
          headerL1   = T.pack . show . IM.size $ m ^. molecule_Atoms
          -- Line 2 of the header contains 1 line of comment. Remove all linebreaks by filtering.
          headerL2   = T.filter (not . isEndOfLine) $ m ^. molecule_Label
          atomLs     =
            IM.foldr' (\atom acc ->
              (T.pack $ printf "%-4s    %12.8F    %12.8F    %12.8F\n"
                (show $ atom ^. atom_Element)
                (fromMaybe 0.0 $ (atom ^. atom_Coordinates) S.!? 0)
                (fromMaybe 0.0 $ (atom ^. atom_Coordinates) S.!? 1)
                (fromMaybe 0.0 $ (atom ^. atom_Coordinates) S.!? 2)
              )
              `T.append`
              acc
            ) "" $ m ^. molecule_Atoms
      in  T.unlines
            [ headerL1
            , headerL2
            ]
          `T.append`
          atomLs

{-|
Write a Tinker XYZ from a 'Molecule'. The writer trusts the '_atom_FFType' to be
correct (if set) and will simply write them out. Therefore it is possible, that wrong atom types can
be written. If they are not set, the writer will simply equalise all atom types to 0, which is OK
for visualisation but obviously not for MM.

This format ingores all deeper level layers of a molecule.
-}
writeTXYZ :: Molecule -> Either String Text
writeTXYZ mol
  | ffTypeCheck =
       toTXYZ <$> (checkMolecule =<< reIndex2BaseMolecule (mol & molecule_SubMol .~ S.empty))
  | otherwise   = Left "writeTXYZ: Not all atoms have Tinker XYZ style atom types."
  where
    -- Check if all atoms have TXYZ atom types.
    ffTypeCheck = all (== TXYZ 0) . IM.map (\a -> a ^. atom_FFType) $ mol ^. molecule_Atoms
    -- Write the molecule as a TXYZ formatted file.
    toTXYZ :: Molecule -> Text
    toTXYZ m =
      let -- The header line contains the number of atoms and separated by a space a comment.
          headerL1 =
            (T.pack . show . IM.size $ mol ^. molecule_Atoms)
            `T.append`
            "  "
            `T.append`
            (T.filter (not . isEndOfLine) $ mol ^. molecule_Label)
            `T.append`
            "\n"
          -- The atom lines contain:
          --   - serial number of the atom
          --   - element symbol of the atom
          --   - XYZ coordinates in angstrom
          --   - The integer number of the atom type
          --   - The indices of all atoms, where to bind to
          atomLs   =
            IM.foldrWithKey' (\key atom acc ->
              -- Print the serial element XYZ FFType
              ( T.pack $ printf "%6d    %2s    %12.8F    %12.8F    %12.8F    %6d"
                  (key + 1)
                  (show $ atom ^. atom_Element)
                  (fromMaybe 0.0 $ (atom ^. atom_Coordinates) S.!? 0)
                  (fromMaybe 0.0 $ (atom ^. atom_Coordinates) S.!? 1)
                  (fromMaybe 0.0 $ (atom ^. atom_Coordinates) S.!? 2)
                  ( case atom ^. atom_FFType of
                      TXYZ i -> i
                      _      -> 0
                  )
              )
              `T.append`
              -- Print the bond targets
              T.stripEnd
              ( IS.foldr' (\x xs ->
                  (T.pack $ printf "    %6d" (x + 1))
                  `T.append`
                  xs
                ) "" $ fromMaybe IS.empty ((m ^. molecule_Bonds) IM.!? key)
              )
              `T.append`
              "\n"
              `T.append`
              acc
            ) "" (m ^. molecule_Atoms)
      in  headerL1
          `T.append`
          atomLs


{-|
Write a simplified .mol2 file (Tripos SYBYL) from a 'Molecule', containing the atoms, connectivities
(single bonds only) and partial charges. This format writes atoms and bonds __only__ from the the
first sublayer of the 'Molecule', which includes the fragment definitions.
-}
writeMOL2 :: Molecule -> Either String Text
writeMOL2 mol
  | ffTypeCheck = toMOL2 <$> (checkMolecule =<< reIndex2BaseMolecule mol)
  | otherwise   = Left "writeMOL2: Not all atoms have MOL2 style atom types."
  where
    ffTypeCheck = all (== Mol2 "") . IM.map (\a -> a ^. atom_FFType) $ mol ^. molecule_Atoms
    toMOL2 :: Molecule -> Text
    toMOL2 m =
      let -- First sublayer of the molecule.
          subMols = m ^. molecule_SubMol
      in  (toMOLECULE m)
           `T.append`
           (toATOM m subMols)
           `T.append`
           (toBOND m)
    -- Write the "@<TRIPOS>MOLECULE" block.
    toMOLECULE :: Molecule -> Text
    toMOLECULE m =
      let nAtoms = IM.size $ m ^. molecule_Atoms
          -- Make the bonds unidirectorial first.
          bonds  = makeBondsUnidirectorial $ m ^. molecule_Bonds
          -- Then get the overall number of bonds.
          nBonds = IM.foldr' (+) 0 . IM.map IS.size $ bonds
      in  T.unlines
            [ "@<TRIPOS>MOLECULE"
            , T.filter (not . isEndOfLine) $ m ^. molecule_Label
            , (T.pack . show $ nAtoms) `T.append` " "
              `T.append`
              (T.pack . show $ nBonds) `T.append` " 0 0 0"
            , "SMALL"
            , "USER_CHARGES"
            , ""
            ]
    -- Write the @<TRIPOS>ATOM block.
    toATOM :: Molecule -> Seq Molecule -> Text
    toATOM m sM =
      let annoFrags = IM.fromAscList . toList . S.zip (S.fromList [ 1 .. ]) $ sM
          atomLines =
            IM.foldrWithKey' (\key atom acc ->
              let -- The index of the submol/fragment, in which the current atom can be found
                  fragmentNum   = findAtomInSubMols key annoFrags :: Maybe Int
                  fragment      = fragmentNum >>= (\fragKey -> IM.lookup fragKey annoFrags)
                  fragmentLabel = _molecule_Label <$> fragment :: Maybe Text
                  thisAtomLine  =
                    T.stripEnd . T.pack $ printf "%7d %-6s %12.8F %12.8F %12.8F %4s %4d  %10s %8s\n"
                      (key + 1)                                           -- Index
                      (show $ atom ^. atom_Label)               -- Atom label
                      (fromMaybe 0.0 $ (atom ^. atom_Coordinates) S.!? 0) -- X
                      (fromMaybe 0.0 $ (atom ^. atom_Coordinates) S.!? 1) -- Y
                      (fromMaybe 0.0 $ (atom ^. atom_Coordinates) S.!? 2) -- Z
                      ( case atom ^. atom_FFType of                       -- Atom label
                          Mol2 t -> T.unpack t
                          _      -> show $ atom ^. atom_Element
                      )
                      (fromMaybe 0 fragmentNum)
                      (fromMaybe "UNL1" $ T.unpack <$> fragmentLabel)
                      (( case atom ^. atom_PCharge of
                          Nothing -> ""
                          Just c  -> printf "%8.4F" c
                      ) :: String)
              in  thisAtomLine `T.append` acc
            ) "" (m ^. molecule_Atoms)
      in  "@<TRIPOS>ATOM\n"
          `T.append`
          atomLines
    -- Write the "@<TRIPOS>BOND" block. All bonds are defined as single bonds for convenience.
    toBOND :: Molecule -> Text
    toBOND m =
      let -- Make the bonds unidirectorial first.
          bonds     = makeBondsUnidirectorial $ m ^. molecule_Bonds
          -- Write a single line for each bond.
          bondLines =
            snd
            $ IM.foldrWithKey' (\origin targets (nthBond, prevBLines) ->
                let -- For each target, write a new line
                    targetLines =
                       snd
                     $ IS.foldr' (\target (nthTarget, prevTLines) ->
                         let targetLine =
                               T.pack $ printf "%6d %6d %6d %4s\n"
                                 (nthTarget + nthBond)
                                 origin
                                 target
                         in  (nthTarget + 1, targetLine `T.append` prevTLines)
                       ) (0, "") targets
                in  (nthBond + IS.size targets, targetLines `T.append` prevBLines)
              ) (1, "") bonds
      in  "@<TRIPOS>BOND" `T.append` bondLines

{-
  "@<TRIPOS>MOLECULE" ++ "\n" ++
  mol ^. molecule_Label ++ "\n" ++
  show nAtoms ++ " " ++ show nBonds ++ " 0 0 0" ++ "\n" ++
  "SMALL" ++ "\n" ++
  "GASTEIGER" ++ "\n\n" ++

  "@<TRIPOS>ATOM" ++ "\n" ++
  concat
    ( map (\(n, a) ->
        printf "%6d    "        n ++
        printf "%-4s    "       (show $ a ^. atom_Element) ++
        (\(x, y, z) -> printf "%12.8F    %12.8F    %12.8F        " x y z)
          (indexAtomCoordinates $ a ^. atom_Coordinates) ++
        printf "%-8s    "
          (if a ^. atom_FFType == ""
             then show (a ^. atom_Element) ++
                  "." ++
                  show (length . I.toList $ a ^. atom_Connectivity)
             else a ^. atom_FFType
          ) ++
        printf "%2d    " (1 :: Int) ++
        printf "%4s    " "UNL1" ++
        printf "%12.8F\n" (fromMaybe 0.0 $ a ^. atom_PCharge)
      ) numberedAtoms
    ) ++ "\n" ++
    "@<TRIPOS>BOND" ++ "\n" ++
    concat
      ( map (\(n, (o, t)) -> printf "    %6d    " n ++
                           printf "%6d    " o ++
                           printf "%6d    " t ++
                           printf "%6d\n" (1 :: Int)
            ) numberedBonds
      )
  where
    atoms = mol ^. molecule_Atoms
    nAtoms = V.length atoms
    numberedAtoms = V.generate nAtoms (\i -> (i, atoms V.! i))
    bonds = map (I.toList . (^. atom_Connectivity)) $ V.toList atoms
    pairBondsRedundant =
      concat
      [ map (\a -> (i + 1, a + 1)) (bonds !! i)
      | i <- [ 0 .. length bonds - 1 ]
      ]
    pairBonds =
      foldr (\a acc -> if swap a `elem` acc
                       then delete a acc
                       else acc
            ) pairBondsRedundant pairBondsRedundant
    nBonds = length pairBonds
    bondIndexList = [ 1 .. nBonds]
    numberedBonds = V.generate nBonds (\i -> (i, pairBonds !! i))
-}

writePDB :: Molecule -> Text
writePDB _mol = undefined

{-|
Write Spicy format, which is an AESON generated JSON document, directly representing the data type
of 'Molecule'.
-}
writeSpicy :: Molecule -> Text
writeSpicy mol = T.decodeUtf8 . encodePretty $ mol
