 -- Convert BlastXML to a directory with one file per query sequence
module BlastCache where

import Data.Binary
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.ByteString.Char8 as BS
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import Data.Vector.Unboxed.Deriving
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable
import System.Directory (createDirectoryIfMissing, removeFile, doesFileExist, doesDirectoryExist)
import Control.DeepSeq
import Control.Arrow (second)
import Data.Hashable
import Text.Printf

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
  -- BlastAlignment::(Bytestring,BlastAlignData)
writeAlignmentCache :: String -> (B.ByteString, [BlastAlignment]) -> IO ()
writeAlignmentCache dir (qsid',ts) = do
  let qsid = B.unpack qsid'
      dnew = printf "%03d" $ hash qsid `mod` 1000
  createDirectoryIfMissing True (dir ++ "/" ++ dnew)
  encodeFile (dir++"/"++dnew++"/"++B.unpack qsid') $ map ( convert.second U.toList) ts
      where convert (name,rest@((A _ _ s):_)) = (name,s,map (\(A a b _) -> (a,b)) rest)


-- | Given the directory and a query sequence, extract the alignments
readAlignmentCache :: String -> String -> IO [BlastAlignment]
readAlignmentCache dir qsid = do
  let dnew = printf "%03d" $ hash qsid `mod` 1000
  dfe <- doesFileExist (dir++"/"++dnew++"/"++qsid)
--  print (dir,dnew,qsid,dfe)
  x <- if dfe then B.readFile (dir++"/"++dnew++"/"++qsid) else B.readFile (dir++"/"++qsid)
--  x <- B.readFile (dir++"/"++qsid)
  let unconvert (name,s,pqs) = let nm = B.copy name
                                   v = U.fromListN (length pqs) $ map (\(p,q) -> (A p q s)) pqs
                               in nm `seq` v `seq` (nm, v)
  return $!! map unconvert $ decode $ x

{-
instance NFData B.ByteString where
  rnf a = a `seq` ()
-}  
{-instance NFData A where
  rnf a = a `seq` ()
-}
