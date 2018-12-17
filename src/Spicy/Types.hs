{-
The module defines all data types that are used in spicy
-}

{-# LANGUAGE DeriveAnyClass  #-}
{-# LANGUAGE DeriveGeneric   #-}
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
, LayerMolecule
, Trajectory
, SuperMolecule
, molecule_Label
, molecule_Atoms
, molecule_Energy
, molecule_Gradient
, molecule_Hessian
--
, Task(..)
, Software(..)
, Basis(..)
, Charge(..)
, Multiplicity(..)
, Shell(..)
, DType(..)
, Efficiency(..)
, has_Gradient
, has_SolventGradient
, has_Hessian
, has_SolventHessian
, SE_Method(..)
, HF_Approx(..)
, DFT_Functional(..)
, DFT_Approx(..)
, MP_Order(..)
, MP_Flavour(..)
, MP_Approx(..)
, CC_Order(..)
, CC_Flavour(..)
, CC_Approx(..)
, CAS_Space(..)
, cas_Space_nElec
, cas_Space_Orbs
, CAS_Roots(..)
, CAS_Weights(..)
, CAS_Approx(..)
--
, QC_SE_Method(..)
, qc_se_Method
, qc_se_Shell
, QC_HF_Shell(..)
, qc_hf_Shell
, qc_hf_Approx
, QC_DFT_Functional(..)
, qc_dft_Functional
, qc_dft_Shell
, qc_dft_Approx
, QC_MPN_Flavour(..)
, qc_mpn_Flavour
, qc_mpn_Shell
, qc_mpn_Approx
, QC_MPN_Order(..)
, qc_mpn_Order
, qc_mpn_Flavour'
, qc_cc_Flavour
, qc_cc_Shell
, qc_cc_Approx
, QC_CC_Order(..)
, qc_cc_Order
, qc_cc_Flavour'
, QC_CAS_Approx(..)
, qc_cas_Approx
, Methods(..)
, methods_qc_SE
, methods_qc_HF
, methods_qc_DFT
, methods_qc_MPN
, methods_qc_CC
, methods_qc_CAS
, dummyMethods
) where
import           Control.DeepSeq
import           Data.IntSet           (IntSet)
import qualified Data.IntSet           as I
import           Data.Map              (Map)
import qualified Data.Map              as Map
import           Data.Maybe
import           Data.Set              (Set)
import qualified Data.Set              as S
import           GHC.Generics          (Generic, Generic1)
import           Lens.Micro.Platform
import           Numeric.LinearAlgebra hiding (Element)
import           Text.Printf


-- | Have a class for the printing of the calculation niveau, which allows to
-- | produce a very readable output
class NiceShow a where
  niceShow :: a -> String    -- printing the isolated object
  niceComplex :: a -> String -- printing the object in the complex form,
                             --   where everything is meant to be printed at once

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
  deriving (Show, Eq, Read, Ord, Enum, Generic, NFData)

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
  , _atom_Connectivity :: IntSet       -- a set of other atoms this one binds to (in the sense of force fields)
                                       --   absolutely meaningless for a single atom, but set on atom level in molecules
                                       --   this is here and not in the molecule layer, because this makes handling with
                                       --   most MM softwares and chemical formats easier (tinker, mol2, PDB)
  } deriving (Eq, Ord, Generic, NFData)
makeLenses ''Atom

-- | A molecule (might be the whole system or a layer, doesnt matter) and all
-- | associated informations
data Molecule = Molecule
  { _molecule_Label    :: String                -- give the molecule a name
  , _molecule_Atoms    :: [Atom]                -- a list of Atoms
  , _molecule_Energy   :: Maybe Double          -- an energy might have been calculated
  , _molecule_Gradient :: Maybe (Vector Double) -- a gradient might have been calculated
  , _molecule_Hessian  :: Maybe (Matrix Double) -- a hessian might have been calculated
  } deriving (Eq, Generic, NFData)
makeLenses ''Molecule

-- | A ONIOM layer with "pseudoatoms" (set 2 atoms in https://doi.org/10.1016/S0166-1280(98)00475-8)
-- | While the molecule is ordinary for this program, pseudo atoms need to be
-- | handled differently from program to program
type LayerMolecule = (Int, Molecule)

-- | Trajectories are simply a list of molecules
type Trajectory = [Molecule]

-- | A supermolecule, which is the whole system (first), and then a list of fragments, treated as
-- | separate molecules. The whole supermolecule containts all atoms and all bonds, but is optional,
-- | as the structure can be completely defined using the fragments. The supermolecule on the other
-- | hand side is suposed to store the results of a calculation (energy, gradient, ...)
type SuperMolecule = (Molecule, [Molecule])


--------------------------------------------------------------------------------
-- Everything that is necessary to describe the calculation niveau
--------------------------------------------------------------------------------
-- | Define commonon tasks we ask the software to perform
data Task =
    Energy
  | Gradient
  | Hessian
  | PartialCharge
  deriving (Show, Eq)
instance NiceShow Task where
  niceShow = show
  niceComplex = show

-- | Software we could talk to
data Software =
    ORCA
  | Psi4
  | NWChem
  | CP2K
  | GAMMES_US
  | Tinker
  deriving (Show, Eq)
instance NiceShow Software where
  niceShow ORCA      = "ORCA"
  niceShow Psi4      = "Psi4"
  niceShow NWChem    = "NWChem"
  niceShow CP2K      = "CP2K"
  niceShow GAMMES_US = "GAMESS US"
  niceShow Tinker    = "Tinker"
  niceComplex = niceShow

-- | Definition of basis sets. Prefer to have a prebuilt list, but allow for the
-- | option to use one, that is not known here by the program specific (!)
-- | keyword (CustomSimple) or by giving the program specific (!) inline notation
-- | of the basis
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
instance NiceShow Basis where
  niceShow Def2_mSVP         = "def2-mSVP"
  niceShow Def2_mTZVP        = "def2-mTZVP"
  niceShow Def2_SV_P         = "def2-SV(P)"
  niceShow Def2_SVP          = "def2-SVP"
  niceShow Def2_TZVP_f       = "def-TZVP(-f)"
  niceShow Def2_TZVP         = "def2-TZVP"
  niceShow Def2_TZVPP        = "def2-TZVPP"
  niceShow Def2_QZVP         = "def2-QZVP"
  niceShow Def2_QZVPP        = "def2-QZVPP"
  niceShow Def2_SVPD         = "def2-SVPD"
  niceShow Def2_TZVPD        = "def2-TZVPD"
  niceShow Def2_TZVPPD       = "def2-TZVPPD"
  niceShow Def2_QZVPD        = "def2-QZVPD"
  niceShow Def2_QZVPPD       = "def2-QZVPPD"
  niceShow Cc_pVDZ           = "cc-pVDZ"
  niceShow Cc_pVTZ           = "cc-pVTZ"
  niceShow Cc_pVQZ           = "cc-pVQZ"
  niceShow Cc_pV5Z           = "cc-pV5Z"
  niceShow Cc_pV6Z           = "cc-pV6Z"
  niceShow Aug_cc_pVDZ       = "aug-cc-pVDZ"
  niceShow Aug_cc_pVTZ       = "aug-cc-pVTZ"
  niceShow Aug_cc_pVQZ       = "aug-cc-pVQZ"
  niceShow Aug_cc_pV5Z       = "aug-cc-pV5Z"
  niceShow Aug_cc_pV6Z       = "aug-cc-pV6Z"
  niceShow (CustomSimple a)  = a
  niceShow (CustomComplex a) = a
  niceComplex = niceShow

-- | Charge and multiplicity of a molecule
type Charge = Int
type Multiplicity = Int

-- | Type of the wavefunction
data Shell =
    Restricted
  | Unrestricted
  | RestrictedOpen
  deriving (Show, Eq)
instance NiceShow Shell where
  niceShow Restricted     = "restricted"
  niceShow Unrestricted   = "unrestricted"
  niceShow RestrictedOpen = "restricted open shell"
  niceComplex Restricted     = "r"
  niceComplex Unrestricted   = "u"
  niceComplex RestrictedOpen = "o"

-- | The way a derivative will be calculated
data DType =
    Analytical
  | Numerical
  | NotAvail
  deriving Eq
instance NiceShow DType where
  niceShow Analytical = "analytical"
  niceShow Numerical  = "numerical"
  niceShow NotAvail   = "not available"
  niceComplex Analytical = "a"
  niceComplex Numerical  = "n"
  niceComplex NotAvail   = "x"

-- | For a given method; how efficient is the calulation of derivatives (beyond
-- | energy)
data Efficiency = Efficiency
  { _has_Gradient          :: [DType]
  , _has_SolventGradient   :: Maybe [DType]
  , _has_Hessian           :: [DType]
  , _has_SolventHessian    :: Maybe [DType]
  , _has_higherDerivatives :: [DType]
  }
makeLenses ''Efficiency
instance NiceShow Efficiency where
  niceShow a =
    "       Gradient        |        Hessian        |     Higher Deriv.     " ++ "\n" ++
    "  GasPhase |  Solvent  |  GasPhase |  Solvent  |                       " ++ "\n" ++
    "-----------+-----------+-----------+-----------+-----------------------" ++ "\n" ++
    printf "  %7s  |  %7s  |  %7s  |  %7s  |  %7s  \n"
      (concatMap niceComplex $ a ^. has_Gradient)
      (case a ^. has_SolventGradient of
        Nothing -> "/"
        Just x  -> concatMap niceComplex x
      )
      (concatMap niceComplex $ a ^. has_Hessian)
      (case a ^. has_SolventHessian of
        Nothing -> "/"
        Just x  -> concatMap niceComplex x
      )
      (concatMap niceComplex $ a ^. has_higherDerivatives)
  niceComplex a =
    "G-" ++
    (concatMap niceComplex $ a ^. has_Gradient) ++
    "(" ++
    (case a ^. has_SolventGradient of
      Nothing -> "/"
      Just x  -> concatMap niceComplex x
    ) ++ ")  " ++
    "H-" ++
    (concatMap niceComplex $ a ^. has_Hessian) ++
    "(" ++
    (case a ^. has_SolventHessian of
      Nothing -> "/"
      Just x  -> concatMap niceComplex x
    ) ++ ")  " ++
    "\n"


  --gradGasPh = a ^. has_Gradient
  --gradSolv = a ^. has_SolventGradient
  --hessGasPh = a ^. has_Hessian
  --hessSolv = a ^. has_SolventHessian


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
  deriving (Eq)
instance Show SE_Method where
  show SE_XTB1 = "GFN-XTB - version 1"
  show SE_XTB2 = "GFN-XTB - version 2"
  show SE_OM1  = "OM(1)"
  show SE_OM2  = "OM(2)"
  show SE_OM3  = "OM(3)"
  show SE_PM3  = "PM3"
  show SE_PM6  = "PM6"
  show SE_PM7  = "PM7"
  show SE_AM1  = "AM1"
  show SE_RM1  = "RM1"
instance NiceShow SE_Method where
  niceShow = show
  --niceComplex = show

-- | Hartree Fock
-- |   -> Approximations
data HF_Approx =
    HF_None
  | HF_RIJK
  | HF_RIJONX
  | HF_RIJCOSX
  deriving (Eq)
instance Show HF_Approx where
  show HF_None    = "None"
  show HF_RIJK    = "RI-JK"
  show HF_RIJONX  = "RI-JONX"
  show HF_RIJCOSX = "RI-JCOSX"

-- | DFT
-- |   -> Functional
-- |   -> Approximations
data DFT_Functional =
  -- LDA
    DFT_PWLDA
  | DFT_VWN3
  | DFT_VWN5
  -- GGA
  | DFT_BLYP
  | DFT_OLYP
  | DFT_GLYP
  | DFT_XLYP
  | DFT_PW91
  | DFT_MPWPW
  | DFT_MPWLYP
  | DFT_PBE
  | DFT_RPBE
  | DFT_REVPBE
  | DFT_PWP
  -- metaGGA
  | DFT_M06_L
  | DFT_M11_L
  | DFT_TPSS
  -- Hybrid
  | DFT_B1LYP
  | DFT_B3LYP
  | DFT_O3LYP
  | DFT_X3LYP
  | DFT_BP86
  | DFT_PBE0
  | DFT_BHANDHLYP
  -- metaHybrid
  | DFT_M06
  | DFT_M06_2X
  | DFT_TPSS0
  | DFT_TPSSh
  -- range separated Hybrid
  | WB97
  | WB97X
  | WB97X_D3
  | CAM_B3LYP
  | LV_BLYP
  -- doubleHybrid
  | DFT_B2PLYP
  | DFT_B2PLYP_RI
  | DFT_B2PLYP_D
  | DFT_B2PLYP_D3
  | DFT_MPW2PLYP
  | DFT_MPW2PLYP_D
  | DFT_B2GP_PLYP
  | DFT_B2K_PLYP
  | DFT_B2T_PLYP
  | DFT_PWPB95
  | DFT_PBPB95_RI
  deriving (Eq, Show)
data DFT_Approx =
    DFT_None
  | DFT_RIJ
  | DFT_RIJK
  | DFT_RIJCOSX
  | DFT_RIJONX
  deriving (Eq)
instance Show DFT_Approx where
  show DFT_None    = "none"
  show DFT_RIJ     = "RI-J"
  show DFT_RIJK    = "RI-JK"
  show DFT_RIJCOSX = "RI-JCOSX"
  show DFT_RIJONX  = "RI-JONX"

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
  deriving (Eq)
instance Show MP_Order where
  show MP_N2   = "2"
  show MP_N2_5 = "2.5"
  show MP_N3   = "3"
  show MP_N4   = "4"
  show MP_N5   = "6"

data MP_Flavour =
    MP_Conv
  | MP_OO
  | MP_F12
  deriving (Eq)
instance Show MP_Flavour where
  show MP_Conv = "conventional"
  show MP_OO   = "orbital optimised"
  show MP_F12  = "explicitly correlated"

data MP_Approx =
    MP_None
  | MP_RI
  | MP_DLPNO
  deriving (Eq)
instance Show MP_Approx where
  show MP_None  = "none"
  show MP_RI    = "RI"
  show MP_DLPNO = "DLPNO"

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
  deriving (Eq)
instance Show CC_Order where
  show CC_D    = "CCD"
  show CC_SD   = "CCSD"
  show CC_SDT  = "CCSDT"
  show CC_SDTQ = "CCSDTQ"
  show CC_SDt  = "CCSD(T)"
  show CC_SDtq = "CCSD(TQ)"

data CC_Flavour =
    CC_Conv
  | CC_OO
  | CC_LR
  | CC_CR
  | CC_F12
  deriving (Eq)
instance Show CC_Flavour where
  show CC_Conv = "conventional"
  show CC_OO   = "orbital optimised"
  show CC_LR   = "locally renormalised"
  show CC_CR   = "completely renormalised"
  show CC_F12  = "explicitly correlated"

data CC_Approx =
    CC_None
  | CC_RI
  | CC_DLPNO
  deriving (Eq)
instance Show CC_Approx where
  show CC_None  = "none"
  show CC_RI    = "RI"
  show CC_DLPNO = "DLPNO"

-- | CAS
-- |   -> Size of the active space. Number of Electrons and List of Orbitals
-- |   -> For state averaging the number of roots
-- |   -> A list of weights of the roots
-- |   -> Approximations applied to the CAS
data CAS_Space = CAS_Space
  { _cas_Space_nElec :: Int
  , _cas_Space_Orbs  :: [Int]
  } deriving (Eq)
makeLenses ''CAS_Space
instance Show CAS_Space where
  show a =
    printf "%10s | %10s\n" "electrons" "orbitals" ++
    printf "%10d | %10d %s"   (a ^. cas_Space_nElec) (length $ a ^. cas_Space_Orbs) (show $ a ^. cas_Space_Orbs)
type CAS_Roots = Int
type CAS_Weights = [Double]
data CAS_Approx =
    CAS_Conv
  | CAS_RIJK
  | CAS_RIJONX
  | CAS_RIJCOSX
  | CAS_DMRG
  deriving (Show, Eq)

-- | A data type, which is designed only for holding quantum chemical
-- | capabilities. It is a highly linked data type containing hierarchical list
-- | of possible QC combination
data QC_SE_Shell = QC_SE_Shell
  { _qc_se_Shell      :: Shell
  , _qc_se_Efficiency :: Efficiency
  }
makeLenses ''QC_SE_Shell
data QC_SE_Method = QC_SE_Method
  { _qc_se_Method :: SE_Method
  , _qc_se_Shell' :: [QC_SE_Shell]
  }
makeLenses ''QC_SE_Method
--------------------------------------------------------------------------------
data QC_HF_Approx = QC_HF_Approx
  { _qc_hf_Approx     :: HF_Approx
  , _qc_hf_Efficiency :: Efficiency
  }
makeLenses ''QC_HF_Approx
data QC_HF_Shell = QC_HF_Shell
  { _qc_hf_Shell   :: Shell
  , _qc_hf_Approx' :: [QC_HF_Approx]
  }
makeLenses ''QC_HF_Shell
--------------------------------------------------------------------------------
data QC_DFT_Approx = QC_DFT_Approx
  { _qc_dft_Approx     :: DFT_Approx
  , _qc_dft_Efficiency :: Efficiency
  }
makeLenses ''QC_DFT_Approx
data QC_DFT_Functional = QC_DFT_Functional
  { _qc_dft_Functional :: DFT_Functional
  , _qc_dft_Shell      :: [Shell]
  , _qc_dft_Approx'    :: [QC_DFT_Approx]
  }
makeLenses ''QC_DFT_Functional
--------------------------------------------------------------------------------
data QC_MPN_Flavour = QC_MPN_Flavour
  { _qc_mpn_Flavour :: MP_Flavour
  , _qc_mpn_Shell   :: [Shell]
  , _qc_mpn_Approx  :: [MP_Approx]
  }
makeLenses ''QC_MPN_Flavour
data QC_MPN_Order = QC_MPN_Order
  { _qc_mpn_Order    :: MP_Order
  , _qc_mpn_Flavour' :: [QC_MPN_Flavour]
  }
makeLenses ''QC_MPN_Order
--------------------------------------------------------------------------------
data QC_CC_Approx = QC_CC_Approx
  { _qc_cc_Approx     :: CC_Approx
  , _qc_cc_Efficiency :: Efficiency
  }
makeLenses ''QC_CC_Approx
data QC_CC_Flavour = QC_CC_Flavour
  { _qc_cc_Flavour :: CC_Flavour
  , _qc_cc_Shell   :: [Shell]
  , _qc_cc_Approx' :: [QC_CC_Approx]
  }
makeLenses ''QC_CC_Flavour
data QC_CC_Order = QC_CC_Order
  { _qc_cc_Order    :: CC_Order
  , _qc_cc_Flavour' :: [QC_CC_Flavour]
  }
makeLenses ''QC_CC_Order
--------------------------------------------------------------------------------
data QC_CAS_Approx = QC_CAS_Approx
  { _qc_cas_Approx     :: CAS_Approx
  , _qc_cas_Efficiency :: Efficiency
  }
makeLenses ''QC_CAS_Approx
--------------------------------------------------------------------------------

data Methods = Methods
  { _methods_qc_SE  :: [QC_SE_Method]
  , _methods_qc_HF  :: [QC_HF_Shell]
  , _methods_qc_DFT :: [QC_DFT_Functional]
  , _methods_qc_MPN :: [QC_MPN_Order]
  , _methods_qc_CC  :: [QC_CC_Order]
  , _methods_qc_CAS :: [QC_CAS_Approx]
  }
makeLenses ''Methods

-- | Generate a string for the possible different wavefunction types
s :: [Shell] -> String
s shell = "(" ++ concatMap shellString shell ++ ") "
  where
    shellString :: Shell -> String
    shellString a
      | Restricted == a = "r"
      | Unrestricted == a = "u"
      | RestrictedOpen == a = "o"
      | otherwise = ""
hfShell :: Shell -> String
hfShell shell
  | shell == Restricted = "RHF"
  | shell == Unrestricted = "UHF"
  | shell == RestrictedOpen = "ROHF"

-- | Generate short strings for ability to use derivatives
d :: Efficiency -> String
d e = gGP ++ gSolv ++ hGP ++ hSolv
  where
    gGP = "G" ++ concatMap dString (e ^. has_Gradient)
    e_gSolv = e ^. has_SolventHessian
    gSolv =
      "(" ++
      case e_gSolv of
        Nothing -> ""
        Just x  -> concatMap dString x
      ++ ") "
    hGP = "H" ++ concatMap dString (e ^. has_Hessian)
    e_hSolv = e ^. has_SolventHessian
    hSolv =
      "(" ++
      case e_hSolv of
        Nothing -> ""
        Just x  -> concatMap dString x
      ++ ") "
    dString a
      | a == Analytical = "a"
      | a == Numerical = "n"
      | a == NotAvail = "/"

dummyMethods = Methods
  { _methods_qc_SE  = []
  , _methods_qc_HF  = []
  , _methods_qc_DFT = []
  , _methods_qc_MPN = []
  , _methods_qc_CC  = []
  , _methods_qc_CAS = []
  }

--------------------------------------------------------------------------------
-- test Types for priting. To be removed later
--------------------------------------------------------------------------------
{-
testSE =
  [ QC_SE_Method
      { _qc_se_Method = SE_PM6
      , _qc_se_Shell' =
          [ QC_SE_Shell
              { _qc_se_Shell = Restricted
              , _qc_se_Efficiency =
                  Efficiency
                    { _has_Gradient = [Analytical]
                    , _has_SolventGradient = Nothing
                    , _has_Hessian = [Numerical]
                    , _has_SolventHessian = Nothing
                    }
              }
          ]
      }
  ] :: [QC_SE_Method]


testHF =
  [ QC_HF_Shell
      { _qc_hf_Shell = RestrictedOpen
      , _qc_hf_Approx =
          [ HF_None
          , HF_RIJK
          , HF_RIJCOSX
          , HF_RIJONX
          ]
      }
  ] :: [QC_HF_Shell]

testDFT =
  [ QC_DFT_Functional
      { _qc_dft_Functional = DFT_TPSS
      , _qc_dft_Shell = [Restricted, Unrestricted, RestrictedOpen]
      , _qc_dft_Approx =
          [ DFT_RIJ
          , DFT_RIJCOSX
          ]
      }
  , QC_DFT_Functional
      { _qc_dft_Functional = DFT_M06
      , _qc_dft_Shell = [Restricted, Unrestricted, RestrictedOpen]
      , _qc_dft_Approx =
          [ DFT_RIJONX
          , DFT_RIJK
          , DFT_RIJCOSX
          ]
      }
  ] :: [QC_DFT_Functional]

testMPN =
  [ QC_MPN_Order
      { _qc_mpn_Order = MP_N2
      , _qc_mpn_Flavour' =
          [ QC_MPN_Flavour
              { _qc_mpn_Flavour = MP_OO
              , _qc_mpn_Shell = [Restricted, Unrestricted]
              , _qc_mpn_Approx =
                  [ MP_None
                  , MP_RI
                  , MP_DLPNO
                  ]
              }
          , QC_MPN_Flavour
              { _qc_mpn_Flavour = MP_Conv
              , _qc_mpn_Shell = [Restricted, Unrestricted]
              , _qc_mpn_Approx =
                  [ MP_RI
                  , MP_DLPNO
                  ]
              }
          ]
      }
  , QC_MPN_Order
      { _qc_mpn_Order = MP_N2_5
      , _qc_mpn_Flavour' =
          [ QC_MPN_Flavour
              { _qc_mpn_Flavour = MP_OO
              , _qc_mpn_Shell = [Restricted, Unrestricted]
              , _qc_mpn_Approx =
                  [ MP_RI
                  , MP_DLPNO
                  ]
              }
          , QC_MPN_Flavour
              { _qc_mpn_Flavour = MP_Conv
              , _qc_mpn_Shell = [Restricted, Unrestricted]
              , _qc_mpn_Approx =
                  [ MP_None
                  ]
              }
          ]
      }
  ]

testCC =
  [ QC_CC_Order
      { _qc_cc_Order = CC_SD
      , _qc_cc_Flavour' =
          [ QC_CC_Flavour
              { _qc_cc_Flavour = CC_Conv
              , _qc_cc_Shell = [Restricted, Unrestricted]
              , _qc_cc_Approx =
                  [ CC_None
                  , CC_RI
                  ]
              }
          ]
      }
  ]

testCAS =
  [ QC_CAS_Approx
      { _qc_cas_Approx = CAS_RIJK }
  ]

testMethods = Methods
  { _methods_qc_SE = testSE
  , _methods_qc_HF = testHF
  , _methods_qc_DFT = testDFT
  , _methods_qc_MPN = testMPN
  , _methods_qc_CC = testCC
  , _methods_qc_CAS = testCAS
  }
-}
