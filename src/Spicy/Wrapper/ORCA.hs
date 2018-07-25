{-
The module defines all data types that are used in spicy
-}
module Spicy.Wrapper.ORCA
(

) where
import           Data.Attoparsec.Text.Lazy
import           Data.Either
import           Data.Map                   (Map)
import qualified Data.Map                   as Map
import           Data.Text                  (Text)
import qualified Data.Text                  as T
import qualified Data.Text.IO               as T
import           Lens.Micro.Platform
import           Spicy.Types
import           Spicy.Wrapper
import           Spicy.Wrapper.ORCA.Methods (methods)
import           System.IO
import           System.Process
