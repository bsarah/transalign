{-# Language DeriveDataTypeable #-}
module Options where

import System.Console.CmdArgs
import Control.Monad (when)
import System.IO (hFlush, hPutStr, stderr)

data Opts = Opts 
            { long, cache :: Bool
            , limit :: Maybe Int
            , blastfiles :: [FilePath]
            , outfile :: Maybe FilePath
            , extract :: [String]
            , cite :: Bool
            , blastfilter :: Maybe Double
            , debug :: Bool
            } deriving (Data,Typeable, Show, Eq)

optdef :: Ann
optdef = opt ("-"::String)

opts :: Opts
opts = Opts
  { long      = False &= help "long output detailing matching positions" &= name "l"
  , limit     = Nothing &= help "max number of alignments to consider" &= name "n"
  , outfile   = Nothing &= help "output to file"
  , extract   = [] &= help "explicit list of alignments to extract"
  , blastfiles = [] &= args &= typ "BLASTXML FILES"
  , cache     = False &= help "generate alignment cache for initial alignment" &= name "c"
  , cite      = False &= help "output citation information"
  , blastfilter = Nothing &= help "exclude intermediate alignment with per-column score less than this value (not using this option disables the filter)"
  , debug     = False &= help "show debug output"
  } 
  &= verbosity
  &= summary "transalign v0.1, ©2012 Ketil Malde"
  &= help "Calculate sequence alignments transitively."
  &= program "transalign"

-- getArgs :: IO Opts
getArgs :: IO ((String -> IO (), String -> IO ()), Opts)
getArgs = do
  o <- cmdArgs opts
  if cite o then error citation else return ()
  loud <- isLoud 
  normal <- isNormal
  let log str = if normal then hPutStr stderr str >> hFlush stderr else return ()
      warn str = if loud then hPutStr stderr str >> hFlush stderr else return ()
  when loud (print o)
  when (length (blastfiles o) /= 2) usage
  return ((warn,log),o)

usage :: t
usage = error "Usage: transalign [options] blast1.xml blast2.xml\n   Use 'transalign --help' for more information."

citation :: String
citation = "\n"++
  "Increasing Sequence Search Sensitivity with Transitive Alignments\n"++
  "Ketil Malde and Tomasz Furmanek\n"++
  "PLoS ONE, 8(2): e54422 (2013)\n"++
  "http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0054422\n"
