{-# Language DoAndIfThenElse  #-}
{-# Language TupleSections    #-}
{-# Language FlexibleContexts #-}

module Main where

import Prelude hiding (log)

import Control.Monad (forM_)
import Control.Monad (when,zipWithM,zipWithM_,forM)
import Control.Parallel.Strategies
import Data.List (sortBy, groupBy, sort, nub)
import Data.Maybe (catMaybes)
import Data.Ord
import GHC.Conc
import GHC.Exts (Down(..))
import GHC.Float (double2Float)
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Map.Strict as M
import qualified Data.Set as S
import qualified Data.Vector.Generic as VG
import System.Directory (doesDirectoryExist,doesFileExist,getDirectoryContents)
import System.Exit
import System.IO.Unsafe
import System.Mem (performGC)
import Text.Printf

import Align (collect_aligns,merge_aligns,A(..))
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
      process_align (debug opts) (blastfilter opts) output sp (B.pack i,myhits)
  else do
    --print up
    myhits <- readAlignments up
    --print $ length myhits
    mapM_ (process_align (debug opts) (blastfilter opts) output sp) myhits

maybeBuildCache :: (String -> IO (),String -> IO ()) -> String -> IO ()
maybeBuildCache (warn,log) sp = do
  d <- doesDirectoryExist (sp++".d")
  if not d
    then do
      log ("Building cache for BLAST output "++show sp++", this might take a while..")
      buildCache sp
      log "..done!\n"
    else do
      l <- doesFileExist (sp++".d"++"/_lock")
      when l $ do
        warn "Lockfile (_lock) in directory indicates caching in progress!\n"
        warn "Either wait, or remove the lock.  Terminating...\n"
        exitWith (ExitFailure 99)

-- process_align :: FilePath -> BlastAlignment -> IO ()
process_align :: Bool -> Maybe Double -> (B.ByteString -> [(Float,B.ByteString,BlastAlignData)] -> IO ()) -> String -> (B.ByteString, [BlastAlignment]) -> IO ()
process_align dbg bfltr output spdir (q, hits) = do
  let hits' = groupBy ((==) `on` fst) . sortBy (comparing fst) $ hits
  let lhits' = length hits'
  let lhits  = length hits
  -- foreach query, collect all target names for streamling hit building
  tgthits <- (M.fromListWith S.union . concat) <$>
    zipWithM (\ kkk hs@((hitname,_):_) -> do
      ts <- readAlignmentCache (spdir++".d") (B.unpack hitname)
      let ttt = unique $ map fst ts
      return (map (,S.singleton hitname) ttt) ) [1 :: Int .. ] hits'
  let ltgthits = M.size tgthits
  -- for each target, run the algorithm
  let tgthitlist = {- take 2 $ -} sortBy (compare `on` (Down . S.size . snd)) $ M.toList tgthits
  tgtgroup <- forM (zip [1 :: Int ..] tgthitlist) $ \(lll, (tgt,hitsources)) -> do
    let lhss = S.size hitsources
    -- foreach query, collect all targets' hits from the cache
    -- let hitnames = unique $ map fst hits
    sphits <- zipWithM ( \ kkk hs@((hitname,_):_) -> unsafeInterleaveIO $ do
                if (S.member hitname hitsources)
                then do
                  when dbg $ printf "# %3d / %3d tgt: %s %5d / %5d [%5d] src: %s\n"
                                    lll ltgthits
                                    (B.unpack tgt)
                                    kkk lhits' lhss
                                    (B.unpack hitname)
                  ts' <- filterAlignments (double2Float <$> bfltr) <$> readAlignmentCache (spdir++".d") (B.unpack hitname)
                  let ts = filter ((==tgt) . fst) ts'
                  let as = {- seq (ts `using` parBuffer numCapabilities rdeepseq) $ -} collect_aligns (const ts) hs
                  return as
                else do
                  return []
              ) [1 :: Int ..] hits'
    return $ merge_aligns $ concat sphits -- (sphits `using` parBuffer numCapabilities rdeepseq)
    
  output q $ reverse
           $ sortBy (compare `on` fst')
           $ map add_score
           $ concat
           $ tgtgroup -- (tgtgroup `using` {- parList rseq) -} parBuffer (2 * numCapabilities) rseq)
--           $ merge_aligns
--           $ concat
--           $ tgtgroup
    where fst' (x,_,_) = x
          f `on` g = \x y -> f (g x) (g y)
  
unique :: Ord a => [a] -> [a]
unique = S.toList . S.fromList

-- | Filter blast alignment to accept only scores @>=f@, where @f@ should
-- be somewhere in the @[0..1]@ range. Values for the per column scores are
-- most likely around @[0.0 .. 0.1]@ but can go higher than @2.0@.
--
-- We assume that each alignment vector only contains entries with the same
-- score values.

filterAlignments :: Maybe Float -> [BlastAlignment] -> [BlastAlignment]
filterAlignments Nothing = id
filterAlignments (Just f) = catMaybes . map go
  where go (s,d)
          | VG.null d = Nothing
          | z>=f      = Just (s,d)
          | otherwise = Nothing
          where A _ _ z = VG.head d -- (s, VG.filter (\(A x y z) -> z >= f) d)

