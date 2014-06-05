import BlastCache (readAlignmentCache)
import Blast (print_alignment)
import Output

import System.Environment (getArgs)
import System.FilePath (splitFileName)
import Data.ByteString.Lazy.Char8 (pack)
import Data.List (partition)

main :: IO ()
main = do 
  (opts,files) <- partition withMinus `fmap` getArgs
  let long = "-l" `elem` opts
  if null files
    then error "No cache files specified"
    else mapM_ (read1 long) files

read1 :: Bool -> FilePath -> IO ()
read1 True fp = do 
  let (d,f) = splitFileName fp
  mapM_ (print_alignment $ pack f) =<< readAlignmentCache d f
read1 False fp = do
  let (d,f) = splitFileName fp
  (output_short $ pack f) =<< map add_score `fmap` readAlignmentCache d f

withMinus :: String -> Bool
withMinus ('-':_) = True
withMinus _ = False
