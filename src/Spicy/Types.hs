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
, atom_PCharge
, atom_Coordinates
, atom_Connectivity
, Molecule(..)
, molecule_Label
, molecule_Atoms
, molecule_Energy
, molecule_Gradient
, molecule_Hessian
--
, Task
, Software
, Basis
, Charge
, Multiplicity
, SE_Method
, HF_Approx
, MP_Order
, MP_Flavour
, MP_Approx
, CC_Order
, CC_Flavour
, CC_Approx
, CAS_Space
, CAS_Roots
, CAS_Weights
, CAS_Approx
) where
import           Data.Map              (Map)
import qualified Data.Map              as Map
import           Lens.Micro.Platform
import           Numeric.LinearAlgebra hiding (Element)

--------------------------------------------------------------------------------
-- Representation of a molecule and its parts
--------------------------------------------------------------------------------
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
  deriving (Show, Eq, Read, Ord, Enum)

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
                                       --   this is here and not in the molecule layer, because this makes handling with
                                       --   most MM softwares and chemical formats easier (tinker, mol2, PDB)
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
-- | While the molecule is ordinary for this program, pseudo atoms need to be
-- | handled differently from program to program
type LayerMolecule = (Int, Molecule)


--------------------------------------------------------------------------------
-- Available calulations implemented in the wrapped software
--------------------------------------------------------------------------------
data Task =
    Energy
  | Gradient
  | Hessian
  | PartialCharge
  deriving (Show, Eq)

data Software =
    ORCA
  | Psi4
  | NWChem
  | Tinker
  | CP2K
  | GAMMES
  deriving (Show, Eq)

data Basis =
  -- Ahlrich
    Def2_mSVP
  | Def2_mTZVP
  | Def2_SV_P
  | Def2_SVP
  | Def2_TZVP_f
  | Def2_TZVP
  | Def2_TZVPP
  | Def2_QZVP
  | Def2_QZVPP
  | Def2_SVPD
  | Def2_TZVPD
  | Def2_TZVPPD
  | Def2_QZVPD
  | Def2_QZVPPD
  -- Dunning
  | Cc_pVDZ
  | Cc_pVTZ
  | Cc_pVQZ
  | Cc_pV5Z
  | Cc_pV6Z
  | Aug_cc_pVDZ
  | Aug_cc_pVTZ
  | Aug_cc_pVQZ
  | Aug_cc_pV5Z
  | Aug_cc_pV6Z
  -- Custom Basis Set. Either as a single keyword or multiline (contraction)
  -- coefficients and so on listed
  | CustomSimple String
  | CustomComplex String
  deriving (Show, Eq)

-- | Charge and multiplicity of a molecule
type Charge = Int
type Multiplicity = Int

-- | Semiempiricism
-- |   -> Available semiempirical hamiltonians
data SE_Method =
    SE_XTB1
  | SE_XTB2
  | SE_OM1
  | SE_OM2
  | SE_OM3
  | SE_PM3
  | SE_PM6
  | SE_PM7
  | SE_AM1
  | SE_RM1
  deriving Eq

-- | Hartree Fock
-- |   -> Approximations
data HF_Approx =
    HF_None
  | HF_RIJK
  | HF_RIJONX
  | HF_RIJCOSX deriving Eq

-- | Møller-Plesset pertubation theory
-- |  -> Order of pertubation theory
-- |  -> Flavour of the MPn. Conventional, orbital optimised or explicitly
-- |     correlated
-- |  -> Approximative MPn calculation
data MP_Order =
    MP_N2
  | MP_N2_5
  | MP_N3
  | MP_N4
  | MP_N5
  deriving Eq
data MP_Flavour =
    MP_Conv
  | MP_OO
  | MP_F12 deriving Eq
data MP_Approx =
    MP_None
  | MP_RI
  | MP_DLPNO
  deriving Eq

-- | Coupled Cluster
-- |  -> Coupled Cluster truncation level
-- |  -> Flavour of the CC. Conventional, orbital optimised, locally
-- |     renormalized, completely renormalized, orbital optimised or explicitly
-- |     correlated
-- |  -> Approximations to the CC.
data CC_Order =
    CC_D
  | CC_SD
  | CC_SDT
  | CC_SDTQ
  | CC_SDt
  | CC_SDtq
  deriving Eq
data CC_Flavour =
    CC_Conv
  | CC_OO
  | CC_LR
  | CC_CR
  | CC_F12
  deriving Eq
data CC_Approx =
    CC_None
  | CC_RI
  | CC_DLPNO
  deriving Eq

-- | CAS
-- |   -> Size of the active space. Number of Electrons and List of Orbitals
-- |   -> For state averaging the number of roots
-- |   -> A list of weights of the roots
-- |   -> Approximations applied to the CAS
type CAS_Space = (Int, [Int])
type CAS_Roots = Int
type CAS_Weights = [Double]
data CAS_Approx =
    CAS_Conv
  | CAS_RIJK
  | CAS_RIJONX
  | CAS_RIJCOSX
  | CAS_DMRG
  deriving Eq
