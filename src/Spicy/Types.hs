{-
The module defines all data types that are used in spicy
-}

{-# LANGUAGE TemplateHaskell #-}
module Spicy.Types
( Element(..)
, AtomLabel
, FFType
, R3Vec
, Atom(..)
, atom_Element
, atom_Label
, atom_IsPseudo
, atom_FFType
, atom_Coordinates
, atom_Connectivity
, Molecule(..)
, molecule_Label
, molecule_Atoms
, molecule_Energy
, molecule_Gradient
, molecule_Hessian
) where
import           Data.Map              (Map)
import qualified Data.Map              as Map
import           Lens.Micro.Platform
import           Numeric.LinearAlgebra hiding (Element)


-- | All chemical elements, have them very clear because force fields and pdb
-- | names may interfer and are just arbitrary strings
data Element =
  H   |                                                                                                                                                                                     He  |
  Li  | Be                                                                                                                                                  | B   | C   | N   | O   | F   | Ne  |
  Na  | Mg                                                                                                                                                  | Al  | Si  | P   | S   | Cl  | Ar  |
  K   | Ca                                                                                      | Sc  | Ti  | V   | Cr  | Mn  | Fe  | Co  | Ni  | Cu  | Zn  | Ga  | Ge  | As  | Se  | Br  | Kr  |
  Rb  | Sr                                                                                      | Y   | Zr  | Nb  | Mo  | Tc  | Ru  | Rh  | Pd  | Ag  | Cd  | In  | Sn  | Sb  | Te  | I   | Xe  |
  Cs  | Ba  | La  | Ce  | Pr  | Nd  | Pm  | Sm  | Eu  | Gd  | Tb  | Dy  | Ho  | Er  | Tm  | Yb  | Lu  | Hf  | Ta  | W   | Re  | Os  | Ir  | Pt  | Au  | Hg  | Tl  | Pb  | Bi  | Po  | At  | Rn  |
  Fr  | Ra  | Ac  | Th  | Pa  | U   | Np  | Pu  | Am  | Cm  | Bk  | Cf  | Es  | Fm  | Md  | No  | Lr  | Rf  | Db  | Sg  | Bh  | Hs  | Mt  | Ds  | Rg  | Cn  | Uut | Fl  | Uup | Lv  | Uus | Uuo
  deriving (Show, Eq, Read, Ord)

-- | An atom label
-- | They may come from pdb or force field parameter files or can be assigned by
-- | other ways just to distinguish specific atoms
type AtomLabel = String

-- | These are labels for molecular mechanics software.
-- | The strings are basically arbitrary and depending on the MM software used
type FFType = String

-- | A vector in the R3 space, used for storing coordinates of atoms (cartesian)
type R3Vec = (Double, Double, Double)

-- | An Atom in a molecule
data Atom = Atom
  { _atom_Element      :: Element      -- element of the atom
  , _atom_Label        :: AtomLabel    -- label, e.g. from a pdb, just for identification, can be empty
  , _atom_IsPseudo     :: Bool         -- boolean telling if this is a pseudo atom, introduced because a bond was broken
  , _atom_FFType       :: FFType       -- label depending on the mm software used, identifying topological atom
  , _atom_PCharge      :: Maybe Double -- possibly a partial charge
  , _atom_Coordinates  :: R3Vec        -- coordinates of the atom, cartesian in R³
  , _atom_Connectivity :: [Int]        -- a list of other atoms this one binds to (in the sense of force fields)
                                       --   absolutely meaningless for a single atom, but set on atom level in molecules
  } deriving (Show, Eq)
makeLenses ''Atom

-- | A molecule (might be the whole system or a layer, doesnt matter) and all
-- | associated informations
data Molecule = Molecule
  { _molecule_Label    :: String                -- give the molecule a name
  , _molecule_Atoms    :: [Atom]                -- a set of Atoms
  , _molecule_Energy   :: Maybe Double          -- an energy might have been calculated
  , _molecule_Gradient :: Maybe (Vector Double) -- a gradient might have been calculated
  , _molecule_Hessian  :: Maybe (Matrix Double) -- a hessian might have been calculated
  } deriving Show
makeLenses ''Molecule

-- | A ONIOM layer with "pseudoatoms" (set 2 atoms in https://doi.org/10.1016/S0166-1280(98)00475-8)
-- | While the molecule is ordinary for this program,
type LayerMolecule = (Int, Molecule)
