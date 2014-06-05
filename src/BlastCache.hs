-- Convert BlastXML to a directory with one file per query sequence
module BlastCache where

import Data.Binary
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.ByteString.Char8 as BS
import System.Directory (createDirectoryIfMissing, removeFile)
import Control.DeepSeq

import Blast (readAlignments, BlastAlignment)
import Align (A(..))

-- | Take a Blast XML file, and split it into a directory with one file
--   per query sequence containing a list of target sequences with alignemnts.
buildCache :: FilePath -> IO ()
buildCache xml = do
  let dir = xml ++ ".d"
  createDirectoryIfMissing True dir
  writeFile (dir++"/_lock") ""
  mapM_ (writeAlignmentCache dir) =<< readAlignments xml
  removeFile (dir++"/_lock")

-- | Write a set of alignments to a file in a specified directory
writeAlignmentCache :: String -> (B.ByteString, [BlastAlignment]) -> IO ()
writeAlignmentCache dir (qsid,ts) = encodeFile (dir++"/"++B.unpack qsid) $ map convert ts
  where convert (name,rest@((A _ _ s):_)) = (name,s,map (\(A a b _) -> (a,b)) rest)

-- | Given the directory and a query sequence, extract the alignments
readAlignmentCache :: String -> String -> IO [BlastAlignment]
readAlignmentCache dir qsid = do
  x <- BS.readFile (dir++"/"++qsid)
  let unconvert (name,s,pqs) = name `seq` (name,map (\(p,q) -> (A p q s)) pqs)
  return $!! map unconvert $ decode $ B.fromChunks [x]

instance NFData B.ByteString where
  rnf a = a `seq` ()
  
instance NFData A where
  rnf a = a `seq` ()