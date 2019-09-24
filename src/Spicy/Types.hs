{-|
Module      : Spicy.Types
Description : Data types used in the Spicy program.
Copyright   : Phillip Seeber, 2019
License     : GPL-3
Maintainer  : phillip.seeber@uni-jena.de
Stability   : experimental
Portability : POSIX, Windows

Spicy.Types contains the definition of all classes and data types, that are used in Spicy. Mainly it
takes care of the description of molecules (structure, topology, potential energy surface, ...).
-}
{-# LANGUAGE DeriveAnyClass    #-}
{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE TemplateHaskell   #-}

module Spicy.Types
  (
  -- * Exceptions
  -- $exceptionTypes
    MolLogicException(..)
  , DataStructureException(..)
  , ParserException(..)
  -- * Performance- and Execution-related Types
  -- $performanceTypes
  , Strat(..)
  , AccVector(..)
  , AccMatrix(..)
  -- * Molecules and Atoms
  -- $moleculeTypes
  , Element(..)
  , AtomLabel
  , FFType(..)
  , Atom(..)
  , atom_Element
  , atom_Label
  , atom_IsPseudo
  , atom_FFType
  , atom_PCharge
  , atom_Coordinates
  , Molecule(..)
  , molecule_Label
  , molecule_Atoms
  , molecule_Bonds
  , molecule_SubMol
  , molecule_Energy
  , molecule_Gradient
  , molecule_Hessian
  , Trajectory
  )
where

import           Control.DeepSeq
import           Control.Exception.Safe

import           Data.Aeson
import           Data.Aeson.Encode.Pretty
import qualified Data.Array.Accelerate         as A
import qualified Data.ByteString.Lazy.UTF8     as BL
import           Data.IntMap.Lazy               ( IntMap )
import           Data.IntSet                    ( IntSet )
import           Data.Sequence                  ( Seq )
import qualified Data.Sequence                 as S
import           Data.Text.Lazy                 ( Text )

import           GHC.Generics                   ( Generic )

import qualified Data.IntMap.Lazy              as IM
import           Lens.Micro.Platform     hiding ( (.=) )
import           Prelude                 hiding ( cycle
                                                , foldl1
                                                , foldr1
                                                , head
                                                , init
                                                , last
                                                , maximum
                                                , minimum
                                                , tail
                                                , take
                                                , takeWhile
                                                , (!!)
                                                )

{-
====================================================================================================
-}
{- $exceptionTypes
These types define composable exceptions for Spicy. The types are specific for a class of problems,
that can occur. These exception types are meant to be used in combination with
'Control.Exception.Safe' as described in
<https://tech.fpcomplete.com/blog/2016/11/exceptions-best-practices-haskell>.
-}
{-|
Exception type for operations on 'Molecule's, which lead to a logical error. This can be caused
because some Spicy assumptions are not met for example.
-}
data MolLogicException = MolLogicException { mlExcFunctionName :: String
                                           , mlExcDescription  :: String
                                           }

instance Show MolLogicException where
  show (MolLogicException f e) = "MoleculeLogicException in function \"" ++ f ++ "\":" ++ e

instance Exception MolLogicException

----------------------------------------------------------------------------------------------------
{-|
Exception type for operations on data structures, which are not meeting necessary criteria for the
operation to perform.
-}
data DataStructureException = DataStructureException { dsExcFunctionName :: String
                                                     , dsExcDescription  :: String
                                                     }

instance Show DataStructureException where
  show (DataStructureException f e) = "DataStructureException in function \"" ++ f ++ "\":" ++ e

instance Exception DataStructureException

----------------------------------------------------------------------------------------------------
{-|
Exception type for textual or binary data, that could not be parsed.
-}
data ParserException = ParserException String

instance Show ParserException where
  show (ParserException e) = "ParserException in parser: \"" ++ e ++ "\""

instance Exception ParserException

{-
====================================================================================================
-}
{- $performanceTypes
These are types related to the type of execution of Spicy. Especially 'AccVector' and 'AccMatrix'
are convenience wrappers around Accelerate data types with some added instances.
-}
{-|
Use serial or parallel processing for large data structures. This helps deciding on a per use base,
if to evaluate in serial or parallel, to avoid nested parallelism. Every method employing parallel
operations by Control.Parallel.Strategies should provide this switch in Spicy.
-}
data Strat = Serial | Parallel
  deriving Eq

----------------------------------------------------------------------------------------------------
{-|
'newtype' wrapper to Accelerate's 'A.Vector's.
-}
newtype AccVector a = AccVector { getAccVector :: A.Vector a
                                }
  deriving ( Generic, Show, Eq )

instance ( ToJSON a, A.Elt a ) => ToJSON (AccVector a) where
  toJSON vec =
    let plainVec        = getAccVector vec
        (A.Z A.:. xDim) = A.arrayShape plainVec
        elements        = A.toList plainVec
    in  object ["shape" .= xDim, "elements" .= elements]

instance ( FromJSON a, A.Elt a ) => FromJSON (AccVector a) where
  parseJSON = withObject "AccVector" $ \vec -> do
    xDim     <- vec .: "shape"
    elements <- vec .: "elements"
    return . AccVector $ A.fromList (A.Z A.:. xDim) elements

----------------------------------------------------------------------------------------------------
{-|
'newtype' wrapper to Accelerate's 'A.Matrix's.
-}
newtype AccMatrix a = AccMatrix { getAccMatrix :: A.Matrix a
                                }
  deriving ( Generic, Show, Eq )

instance ( ToJSON a, A.Elt a ) => ToJSON (AccMatrix a) where
  toJSON mat =
    let plainMat                  = getAccMatrix mat
        (A.Z A.:. xDim A.:. yDim) = A.arrayShape plainMat
        elements                  = A.toList plainMat
    in  object ["shape" .= (xDim, yDim), "elements" .= elements]

instance ( FromJSON a, A.Elt a ) => FromJSON (AccMatrix a) where
  parseJSON = withObject "AccVector" $ \mat -> do
    (xDim, yDim) <- mat .: "shape"
    elements     <- mat .: "elements"
    return . AccMatrix $ A.fromList (A.Z A.:. xDim A.:. yDim) elements

{-
====================================================================================================
-}
{- $moleculeTypes
Types to describe a molecule from a computational geometry point of view. The representation of a
molecule aims to be as descriptive as unambigous as possible and capture all quantities relevant for
molecular dynamics and geometry optimisations (energies, gradients, hessians, force fiekd types,
coordinates, connectivities). Electronic properties, such as dipoles, orbitals, the wavefunction and
so on, will be ignored, as Spicy is not aiming to be an analysis program, but a interface to create
multilayer calculations. Analysis of the results is up to postprocessing.
-}
{-|
All chemical elements. Have them very clear because force fields and pdb names may interfer and are
just arbitrary strings.
-}
data Element
  = H
  | He
  | Li
  | Be
  | B
  | C
  | N
  | O
  | F
  | Ne
  | Na
  | Mg
  | Al
  | Si
  | P
  | S
  | Cl
  | Ar
  | K
  | Ca
  | Sc
  | Ti
  | V
  | Cr
  | Mn
  | Fe
  | Co
  | Ni
  | Cu
  | Zn
  | Ga
  | Ge
  | As
  | Se
  | Br
  | Kr
  | Rb
  | Sr
  | Y
  | Zr
  | Nb
  | Mo
  | Tc
  | Ru
  | Rh
  | Pd
  | Ag
  | Cd
  | In
  | Sn
  | Sb
  | Te
  | I
  | Xe
  | Cs
  | Ba
  | La
  | Ce
  | Pr
  | Nd
  | Pm
  | Sm
  | Eu
  | Gd
  | Tb
  | Dy
  | Ho
  | Er
  | Tm
  | Yb
  | Lu
  | Hf
  | Ta
  | W
  | Re
  | Os
  | Ir
  | Pt
  | Au
  | Hg
  | Tl
  | Pb
  | Bi
  | Po
  | At
  | Rn
  | Fr
  | Ra
  | Ac
  | Th
  | Pa
  | U
  | Np
  | Pu
  | Am
  | Cm
  | Bk
  | Cf
  | Es
  | Fm
  | Md
  | No
  | Lr
  | Rf
  | Db
  | Sg
  | Bh
  | Hs
  | Mt
  | Ds
  | Rg
  | Cn
  | Uut
  | Fl
  | Uup
  | Lv
  | Uus
  | Uuo
  deriving ( Show, Eq, Read, Ord, Enum, Generic, NFData )

instance ToJSON Element where
  toEncoding = genericToEncoding defaultOptions

instance FromJSON Element

----------------------------------------------------------------------------------------------------
{-|
An atom label. They may come from pdb or force field parameter files or can be assigned by other
ways just to distinguish specific atoms.
-}
type AtomLabel = Text

----------------------------------------------------------------------------------------------------
{-|
These are labels for molecular mechanics software. The strings are basically arbitrary and depending
on the MM software used.
-}
data FFType = Mol2 Text | TXYZ Int | PDB Text | XYZ
  deriving ( Generic )

makeLenses ''FFType

instance Eq FFType where
  Mol2 _ == Mol2 _ = True
  Mol2 _ == _      = False
  TXYZ _ == TXYZ _ = True
  TXYZ _ == _      = False
  PDB  _ == PDB _  = True
  PDB  _ == _      = False
  XYZ    == XYZ    = True
  XYZ    == _      = False

instance ToJSON FFType where
  toEncoding = genericToEncoding defaultOptions

instance FromJSON FFType

----------------------------------------------------------------------------------------------------
{-|
An Atom in a 'Molecule'. Atoms are compared by their indices only and they must therefore be unique.
The coordinates of the 'Atom' are defined as 'Seq', as this is extremely easy to concatenate when
building a coordinate vector.
-}
data Atom =
  Atom { _atom_Element     :: Element      -- ^ Chemical 'Element' of the atom.
       , _atom_Label       :: AtomLabel    -- ^ Label, e.g. from a pdb, just for identification, can
                                           --   be empty.
       , _atom_IsPseudo    :: Bool         -- ^ Boolean, telling if this is a pseudo atom,
                                           --   introduced because a bond was broken.
       , _atom_FFType      :: FFType       -- ^ Label depending on the MM software used, identifying
                                           --   topological atom.
       , _atom_PCharge     :: Maybe Double -- ^ Possibly a partial charge.
       , _atom_Coordinates :: Seq Double   -- ^ Coordinates of the atom, cartesian in R³. Relies on
                                           --   the parser to fill with exactly 3 values.
       }
  deriving ( Eq, Generic )

makeLenses ''Atom

instance Show Atom where
  show = BL.toString . encodePretty

instance ToJSON Atom where
  toEncoding = genericToEncoding defaultOptions

instance FromJSON Atom

----------------------------------------------------------------------------------------------------
{-|
A molecule, which might be the whole system, an ONIOM layer or a fragment of the system, each
containing possibly even higher layers for ONIOM or fragments. Stores all associated informations of
a layer.

Starting from a top level molecule, all atoms and bonds of the system are expected to be in the in
this top layer (except pseudoatoms of deeper layers). Therefore if atoms are in a deeper layers of
the recursion, their information is not used to describe a higher layer. Instead, all atoms of
deeper layers (except pseudoatoms) must be also replicated in a higher layer.
-}
data Molecule =
  Molecule { _molecule_Label    :: Text                     -- ^ Comment or identifier of a
                                                            --   molecule. Can be empty.
           , _molecule_Atoms    :: IntMap Atom              -- ^ An 'IntMap' of 'Atom's, 'Atom's
                                                            --   identified by their 'Int' index.
           , _molecule_Bonds    :: IntMap IntSet            -- ^ An IntMap, mapping the index of an
                                                            --   'Atom' in the 'Molecule' to the
                                                            --   indices of all 'Atom's, to which
                                                            --   it binds.
           , _molecule_SubMol   :: Seq Molecule             -- ^ A Molecule might contain other
                                                            --   molecules. These might be fragments
                                                            --   or higher level ONIOM layers.
           , _molecule_Energy   :: Maybe Double             -- ^ An energy, that might have been
                                                            --   calculated.
           , _molecule_Gradient :: Maybe (AccVector Double) -- ^ A gradient, that might have been
                                                            --   calculated.
           , _molecule_Hessian  :: Maybe (AccMatrix Double) -- ^ A hessian, that might have been
                                                            --   calculated.
           }
  deriving ( Eq, Generic )

makeLenses ''Molecule

instance Show Molecule where
  show = BL.toString . encodePretty

instance ToJSON Molecule where
  toEncoding = genericToEncoding defaultOptions

instance FromJSON Molecule

----------------------------------------------------------------------------------------------------
{-|
Trajectories are simply 'Seq'uences of 'Molecule's.
-}
type Trajectory = Seq Molecule

{-|
Dummy molecule for testing purposes
-}
dummyMolecule :: Molecule
dummyMolecule = Molecule "" IM.empty IM.empty S.empty Nothing Nothing Nothing
