{-|
Module      : Spicy.Internal.Accelerate
Description : Wrapper for device specific Accelerate/LLVM operations.
Copyright   : Phillip Seeber, 2019
License     : GPL-3
Maintainer  : phillip.seeber@uni-jena.de
Stability   : experimental
Portability : POSIX, Windows

A module to deal with device specific Accelerate/LLVM functionalities on a cabal level. This avoids
using CPP directives, which potentially confuse compilers and prettifiers.
-}
module Spicy.Internal.Accelerate
  ( runQ
  , runN
  )
where
import           Data.Array.Accelerate.LLVM.Native
