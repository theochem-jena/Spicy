module Spicy.Wrapper.ORCA
(

) where
import           Data.Attoparsec.Text.Lazy
import           Data.Either
import           Data.Map                  (Map)
import qualified Data.Map                  as Map
import           Data.Text                 (Text)
import qualified Data.Text                 as T
import qualified Data.Text.IO              as T
--import           Lens.Micro.Platform
import           Spicy.Types
import           Spicy.Wrapper
import           System.Process
import System.IO

runORCA_IO :: FilePath -> IO Text
runORCA_IO inputFile = do
  let tmpORCA = "/tmp"
  -- ask the system for the full orca path with a "which orca"
  (_, Just outPath, Just errPath, _) <-
    createProcess $ (proc "which" ["orca"])
    { std_out = CreatePipe
    , std_err = CreatePipe
    }
  pathORCA <- hGetContents outPath

  (_, Just outORCA, Just errORCA, _) <-
   createProcess $ (proc pathORCA [inputFile])
   { std_out = CreatePipe
   , cwd = Just tmpORCA
   }
  orcaOutput <- T.hGetContents outORCA
  return orcaOutput
