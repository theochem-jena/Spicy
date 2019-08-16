{-|
Module      : Spicy.MolecularSystem
Description : Handling molecular informations
Copyright   : Phillip Seeber, 2019
License     : GPL-3
Maintainer  : phillip.seeber@uni-jena.de
Stability   : experimental
Portability : POSIX, Windows

This module deals with the partitioning of the system, creation of bonds, assignment of
substructeres to layers and creation of ghost atoms. Following conventions shall apply:

    * highest level region has highest index

    * lowest level region has index 0 and contains the complete system
-}
module Spicy.MolecularSystem
(
) where
