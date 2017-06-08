{-# Language DoAndIfThenElse #-}

module Main where

import Prelude hiding (log)

import qualified Data.Map.Strict as M
import qualified Data.Set as S
import System.Directory (doesDirectoryExist,doesFileExist,getDirectoryContents)
import Control.Monad (when)
import System.Exit
import qualified Data.ByteString.Lazy.Char8 as B
import Data.List (sortBy, groupBy)
import Control.Monad (forM_)
import Data.Ord
import System.Mem (performGC)
import Control.Parallel.Strategies

import Align (collect_aligns,merge_aligns)
import Blast (BlastAlignment, BlastAlignData, readAlignments)
import BlastCache (readAlignmentCache, buildCache)
import Options (getArgs, Opts(..))
import Output (output_long, output_short, add_score)

main :: IO ()
main = do
  (log, opts) <- getArgs
  let output = if long opts then output_long
               else output_short
  -- read all alignments vs UniProt
  let [up,sp] = blastfiles opts
  maybeBuildCache log sp
  
  if cache opts then do
    maybeBuildCache log up
    -- uphits <- readAlignments up -- :: [ (query, [(target,alignment)] ) ]
    inputs <- if (not $ null $ extract opts) 
              then return (extract opts) 
              else filter (\x -> x /= "." && x /= "..") `fmap` getDirectoryContents (up++".d")
  
    flip mapM_ inputs $ \i -> do  -- forM_
      myhits <- readAlignmentCache (up++".d") i
      process_align output sp (B.pack i,myhits)
  else do
    --print up
    myhits <- readAlignments up
    --print $ length myhits
    mapM_ (process_align output sp) myhits

maybeBuildCache :: (String -> IO (),String -> IO ()) -> String -> IO ()
maybeBuildCache (warn,log) sp = do
  d <- doesDirectoryExist (sp++".d")
  if not d
    then do
      log ("Building cache for BLAST output "++show sp++", this might take a while..\n")
      buildCache sp
      log "..done!\n"
    else do
      l <- doesFileExist (sp++".d"++"/_lock")
      when l $ do
        warn "Lockfile (_lock) in directory indicates caching in progress!\n"
        warn "Either wait, or remove the lock.  Terminating...\n"
        exitWith (ExitFailure 99)

-- process_align :: FilePath -> BlastAlignment -> IO ()
process_align :: (B.ByteString -> [(Float,B.ByteString,BlastAlignData)] -> IO ()) -> String -> (B.ByteString, [BlastAlignment]) -> IO ()
process_align output spdir (q, hits) = do
  let hits' = groupBy ((==) `on` fst) . sortBy (comparing fst) $ hits
  --print $ length hits
  --print $ length hits'
  -- foreach query, collect all targets' hits from the cache
  -- let hitnames = unique $ map fst hits
  sphits <- mapM ( \ hs@((hitname,_):_) -> do
              --print "before gc"
              --print $ length hs
              performGC
              --print "after gc"
              ts <- readAlignmentCache (spdir++".d") (B.unpack hitname)
              --print "after rac"
              let as = parMap rdeepseq id $ collect_aligns (const ts) hs
              --print "after collect"
              return as) hits'
    
  --  sphits <- M.fromList `fmap` mapM (\h -> do 
  --                                     ts <- readAlignmentCache (spdir++".d") (B.unpack h)
  --                                     return (h,ts)) [hitname]
--  forM_ (M.assocs sphits) $ \ (h,ts) -> do
--    print (h,length ts)
--  let as = merge_aligns $ collect_aligns (mlu sphits) hits
--      mlu m k = maybe [] id $ M.lookup k m
  output q $ reverse $ sortBy (compare `on` fst') $ map add_score $ merge_aligns $ concat sphits 
    where fst' (x,_,_) = x
          f `on` g = \x y -> f (g x) (g y)
  
unique :: Ord a => [a] -> [a]
unique = S.toList . S.fromList

