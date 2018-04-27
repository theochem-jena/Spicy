module Spicy.Wrapper
( Task
, Software
) where

data Task = Energy | Gradient | Hessian | PartialCharge deriving (Show, Eq)

data Software =
  ORCA
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

type Charge = Int
type Multiplicity = Int

-- | Semiempiricism
-- |   -> Available semiempirical hamiltonians
data SEMethod = XTB1 | XTB2 | OM1 | OM2 | OM3 | PM3 | PM6 | PM7 | AM1 | RM1
  deriving Eq

-- | Hartree Fock
-- |   -> Approximations
data HFApprox = None | RIJK | RIJONX | RIJCOSX deriving Eq

-- | Møller-Plesset pertubation theory
-- |  -> Order of pertubation theory
-- |  -> Flavour of the MPn. Conventional, orbital optimised or explicitly
-- |     correlated
-- |  -> Approximative MPn calculation
data MPOrder = N2 | N2_5 | N3 | N4 | N5 deriving Eq
data MPFlavour = Conv | OO | F12 deriving Eq
data MPApprox = None | RI | DLPNO deriving Eq

-- | Coupled Cluster
-- |  -> Coupled Cluster truncation level
-- |  -> Flavour of the CC. Conventional, orbital optimised, locally
-- |     renormalized, completely renormalized, orbital optimised or explicitly
-- |     correlated
-- |  -> Approximations to the CC.
data CCOrder = D | SD | SDT | SDTQ | SD_T | SD_TQ deriving Eq
data CCFlavour = Conv | OO | LR | CR | OO | F12 deriving Eq
data CCApprox = None | RI | DLPNO deriving Eq

-- | CAS
-- |   -> Size of the active space. Number of Electrons and List of Orbitals
-- |   -> For state averaging the number of roots
-- |   -> A list of weights of the roots
-- |   -> Approximations applied to the CAS
type CASSpace = (Int, [Int])
type CASRoots = Int
type CASWeights = [Double]
data CASApprox = Conv | RIJK | RIJONX | RIJCOSX | DMRG deriving Eq



data Method =
    MM Forcefield -- Molecular Mechanics with a specific forcefield
  | SE SEMethod  -- Semiempiricism
  | HF -- Hartree-Fock
  |
  | LCCD
  | CCD -- Coupled Cluster
  | CCSD
  | CCSDT
  | CCSDTQ
  | CCSD_T
