-- | Parser test script
--runghc -package-db --ghc-arg=.cabal-sandbox/x86_64-linux-ghc-8.0.1-packages.conf.d/ BlastParserTest.hs blasttest/BlastPtest1.csv
--   read from file and directly print parsing output

module Main where
    
import System.Environment (getArgs)
import System.Process 
import Text.ParserCombinators.Parsec
import System.IO
import System.Environment
import TabularBlastParser
import Data.Either.Unwrap
import qualified Data.ByteString.Char8 as C
    
main = do
  args <- getArgs
  let inputFile = (head args)
  inputData <- readTabularBlasts inputFile
  let rightBlast = fromRight inputData
  print rightBlast

