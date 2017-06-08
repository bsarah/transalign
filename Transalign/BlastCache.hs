{-# Language OverloadedStrings #-}
{-# Language DoAndIfThenElse #-}

-- Convert BlastXML to a directory with one file per query sequence
module Transalign.BlastCache where

import Data.Binary
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.ByteString.Char8 as BS
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import Data.Vector.Unboxed.Deriving
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable
import System.Directory (createDirectoryIfMissing, removeFile, doesFileExist, doesDirectoryExist)
import System.FilePath ((</>))
import Control.DeepSeq
import Control.Arrow (second)
import Data.Hashable
import Text.Printf
import Debug.Trace
import System.IO (stderr,hPutStrLn)
import qualified Codec.Compression.LZ4 as LZ4
import qualified Control.Concurrent.ParallelIO.Local as PIO
import qualified Data.List.Split as LS
import System.Mem
import qualified Data.Serialize as Cereal

import Transalign.Blast (readAlignments, BlastAlignment)
import Transalign.Align (A(..))

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
--writeAlignmentCache dir (qsid',[]) = return
writeAlignmentCache dir (qsid',ts) = do
  let qsid = B.unpack qsid'
      dnew = printf "%03d" $ hash qsid `mod` 1000
  createDirectoryIfMissing True (dir ++ "/" ++ dnew)
  BS.writeFile (dir++"/"++dnew++"/"++B.unpack qsid')
    $ maybe (error "writeAlignmentCache: LZ4.compress failed!") id
    $ LZ4.compress
    $ Cereal.encode
    $ map ( convert.second U.toList) ts
      where convert (name,rest@((A _ _ s):_)) = (name,s,map (\(A a b _) -> (a,b)) rest)


-- | Given the directory and a query sequence, extract the alignments
readAlignmentCache :: String -> String -> IO [BlastAlignment]
readAlignmentCache dir qsid = do
  let dnew = printf "%03d" $ hash qsid `mod` 1000
  let fnew = dir </> dnew </> qsid
  let fold = dir </> qsid
  dfeNew <- doesFileExist fnew
  dfeOld <- doesFileExist fold
  xNew <- if dfeNew
            then maybe (error "readAlignmentCache: LZ4.decompress failed!") id
                  <$> LZ4.decompress
                  <$> BS.readFile fnew
            else return BS.empty
  xOld <- if dfeOld
            then BS.readFile fold
            else return BS.empty
  let unconvert (name,s,pqs) = let nm = B.fromStrict $ BS.copy name
                                   v = U.fromListN (length pqs) $ map (\(p,q) -> (A p q s)) pqs
                               in nm `seq` v `seq` (nm, v)
  case (BS.null xNew , BS.null xOld) of
    (True,True) -> do hPutStrLn stderr $ "ERROR: cache file not found: " ++ fnew
                      return []
    (True,_   ) -> return $!! map unconvert $ either error id $ Cereal.decode $ xOld
    (_   ,True) -> return $!! map unconvert $ either error id $ Cereal.decode $ xNew
    (_   ,_   ) -> return []

{-
instance NFData B.ByteString where
  rnf a = a `seq` ()
-}  
{-instance NFData A where
  rnf a = a `seq` ()
-}
