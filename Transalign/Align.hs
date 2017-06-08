{- IDEA:  annotate a transcript by doing blastx vs NR, then blastp NR vs SwissProt
   Align transcript vs SP by mapping each position via NR hits to positions in SP.
   How to test?

   Generalize: 
     annotate a_vs_b.xml b_vs_c.xml c_vs_d.xml => result: d_vs_a (or: a_vs_d?)
-}

{-# Language TypeFamilies #-}
{-# Language MultiParamTypeClasses #-}
{-# Language TemplateHaskell #-}
{-# Language FlexibleInstances #-}
-- {-# Language UndecidableInstances #-}
{-# Language FlexibleContexts #-}
{-# LANGUAGE CPP #-}

module Transalign.Align (A(..),Alignment,score,collect_aligns,merge_aligns) where

import Control.Applicative ( (<$>) )
import qualified Data.IntMap.Strict as M
import Data.Hashable
import qualified Data.HashMap.Strict as H
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import Data.Vector.Unboxed.Deriving
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable
import qualified Data.Judy as J
import Data.List (sort,sortBy,foldl',last,groupBy)
import Data.Int (Int32)
import Data.Word (Word32)
import Debug.Trace
import Control.DeepSeq
import Control.Parallel.Strategies
import Data.Vector.Algorithms.Merge
import Control.Monad
import System.IO.Unsafe
import Unsafe.Coerce
import Data.Maybe
import GHC.Conc.Sync
import System.Mem (performGC)
import Control.Seq as CS


type Pos = Int32
type Score = Float

data A  = A {-# UNPACK #-} !Pos  {-# UNPACK #-} !Pos  {-# UNPACK #-} !Score deriving (Show,Eq)

derivingUnbox "A"
    [t| A -> (Pos,Pos,Score) |]
    [| \(A p q s) -> (p,q,s) |]
    [| \(p,q,s) -> A p q s |]


-- A proper 'Alignment' links positions in one sequence to positions in another,
-- assigning a score to each link.  The position links are co-linear and one-to-one.


type Alignment score pos1 pos2 = U.Vector A


-- | get sequences from target which fit to query sequence 
collect_aligns :: (usid -> [(ssid, Alignment s q r)]) -> [(usid, Alignment s p q)] -> [(ssid,Alignment s p r)]
collect_aligns sp_lu ups = do -- list monad
  (uhit,ualign) <- ups
  (shit,salign) <- sp_lu uhit
  return (shit,trans_align ualign salign)


score :: Alignment s p q -> Score -- s
score = G.sum . G.map (\(A _ _ s) -> s)

-- | Reverse the alignment
invert :: Alignment s p q -> Alignment s q p
invert = G.map (\(A x y s) -> (A y x s))


-- | Follow positions via the first 'Alignment', through the second to 
--   form the transitive alignment.  Note that if the alignments have 
--   nothing in common, this may produce an empty alignment.
trans_align :: Alignment s p q -> Alignment s q r -> Alignment s p r
trans_align xorig yorig = G.force $ G.unfoldr go (xorig,yorig) where
    go (x,y) 
        | G.null x || G.null y = Nothing
        | otherwise = let (A xp xq xs) = G.unsafeHead x
                          (A yq yr ys) = G.unsafeHead y
                          xt = G.unsafeTail x
                          yt = G.unsafeTail y
                      in
                        if (xq == yq) then Just (A xp yr (min xs ys),(xt,yt))
                        else
                            if (xq < yq) then go(xt,y)
                            else go(x,yt)
{-# NOINLINE trans_align #-}

-- | collect hits against the same sequence, and calculate a consensus alignment
merge_aligns :: (Show ssid, Ord ssid, NFData ssid) => [(ssid,Alignment s p q)] -> [(ssid,Alignment s p q)]
merge_aligns xs =  let y = groups $ filter is_not_empty xs
                   in map merge_group y `CS.using` CS.seqList CS.rdeepseq
    where is_not_empty = not . G.null . snd 


-- | Use group_al, merge1 and sort_column to get a list/vector of possible alignments including
-- overall scores of each alignment. Group_al computes a matrix with (p,[(q,s)]). Now compute the best alignment calculating its overall score
-- the list has to be kept sorted increasingly by scores with sort_column
merge_group :: (Show ssid) => (ssid, V.Vector (Alignment s p q)) -> (ssid, Alignment s p q)
merge_group (tgt,as) = let bs = group_al5 as 
                       in (tgt , go [] bs) -- $ map (\(p,xs) -> (p,collect xs)) 
    where        
    -- this is just a regular alignment using the 
    -- scores from the provided alignments. Sparse matrix.
      go :: [(Score,[A])] -> V.Vector (Pos, U.Vector (Pos,Score)) -> Alignment s p q
      go y x
          | null y && G.null x = error ("last will fail! " ++ show tgt ++ " " ++ show as)
          | null y = let (p,xs) = G.unsafeHead x
                         rest = G.unsafeTail x
                     in go (map (\(q,s)->(s,[(A p q s)])) (U.toList xs)) rest                           
          | G.null x = U.fromList $ reverse $ snd $ last $ y
          | otherwise = let xh = G.unsafeHead x
                            (p,xv) = xh
                            xt = G.unsafeTail x
                        in go (sort_column $ merge1 y xh) xt


    

-- | DP part. p x q matrix with scores as entries. Compute the best overall alignment.
--   x is the 'prev' column, the current results are in prev and compared in each DP-step. 
merge1 :: [(Score, [A])] -> (Pos, U.Vector (Pos, Score)) -> [(Score, [A])] --dynamic programming, forward tracking, not many costs: 
merge1 x y
  | null x && (U.null $ snd y) = error "error in merge1, input 0"
  | null x = error "error in merge1, empty prev"
  | U.null $ snd y = x
  | otherwise = let (ts,as) = head x 
                    ht = tail x 
                    (A p1 q1 s1) = head as
                    ast = tail as 
                    (p,ys) = y
                    (q,s) = G.unsafeHead ys
                    nc = G.unsafeTail ys                  
                in
                  if q <= q1 then (s,[(A p q s)]) : merge1 x (p,nc)
                  else if null ht then (ts+s,(A p q s):(A p1 q1 s1):ast) : merge1 x (p,nc)
                       else let (_,hpc) = head ht
                                (A _ q2 _) = head hpc
                            in if (q2 < q)
                               then merge1 ht y      
                               else (ts+s,(A p q s):(A p1 q1 s1):ast) : merge1 x (p,nc)
  

-- | sort on q-coordinate and filter so scores are increasing
sort_column :: [(Score, [A])] -> [(Score, [A])]
sort_column = filter_scores . Data.List.sortBy (compare `on` qval)
    where qval (sc,as) = let (A p q s) = head as
                         in q
          f `on` g = \x y -> f (g x) (g y)
          filter_scores ((xs,xss):(ys,yss):rest)
            | ys < xs = filter_scores ((xs,xss):rest)
            | otherwise = (xs,xss):filter_scores ((ys,yss):rest)
          filter_scores [x] = [x]
          filter_scores []  = []


        
-- | Add scores for all (q,s) pairs mapping to the same q. Input must be sorted on q.
collect :: [(Pos, Score)] -> [(Pos, Score)]
collect ((q1,s1):(q2,s2):rest)  
  | q1==q2    = collect ((q1,s1+s2):rest)
  | q1 < q2   = (q1,s1) : collect ((q2,s2):rest)
  | otherwise = error ("collect: unsorted input!\n"++show ((q1,s1):(q2,s2):take 2 rest))
collect [x] = [x]
collect []  = []


-- | Group alignments against the same target
-- Basically a quicksort
groups :: Ord sid => [(sid, Alignment s p q)] -> [(sid, V.Vector (Alignment s p q))]
groups ((s,a):xs) = this : groups less ++  groups more
    where this = (s, V.fromList $ a : map snd (filter ((s==).fst) xs))
          less = filter ((s<).fst) xs
          more = filter ((s>).fst) xs
groups [] = []
-- prop: length = length . concat . groups

-- Group alignments on first coordinate (p), using a direct merge of
-- the alignments. (This is "al2" in tests).
-- todo: simultaneously collect in a Map q s?
-- This is slower than group_al'
{-
group_al :: V.Vector (Alignment s p q) -> V.Vector (Pos, V.Vector (Pos,Score))
group_al = join . (V.foldr1 merge)
  where        
    merge one two
      | V.null one = two
      | V.null two = one
      | otherwise = let (A p1 q1 s1) = V.head one
                        one_rest = V.tail one
                        (A p2 q2 s2) = V.head two
                        two_rest = V.tail two
                    in if p1 < p2 then (A p1 q1 s1) `V.cons` merge one_rest two
                       else if p1 == p2 then (A p1 q1 s1) `V.cons` ((A p2 q2 s2) `V.cons` merge one_rest two_rest)
                            else (A p2 q2 s2) `V.cons` merge one two_rest
    join x = let (A p q s) = V.head x
                 xs = V.tail x
             in go p (V.singleton (q,s)) xs
    go p0 acc x 
      | V.null x = V.singleton (p0, sort' acc)
      | otherwise = let (A p q s) = V.head x
                        xs = V.tail x
                    in if p0 == p then go p0 ((q,s) `V.cons` acc) xs 
                       else (p0, sort' acc) `V.cons` join x
    sort' = (V.fromList) . sort . (V.toList)

-}
      {-  sort' :: V.Vector (Pos,Score) -> V.Vector (Pos,Score) -}


{-
-- Using a Map is slightly better than using quicksort OR mergesort(!)
group_al' :: (Ord s) => [Alignment s Int Int] -> [(Int,[(Int,s)])]
group_al' =  merge . concatMap (map (\(p,q,s) -> (p,[(q,s)])))
   where merge = M.toList . M.map sort . M.fromListWith (flip (++))

-- As alignmnents are all sorted already, quicksort (from groups) is a 
-- *very* bad choice.  This is hopelessly inefficient.
group_al'' :: (Ord p, Ord q, Ord s) => [Alignment s p q] -> [(p,[(q,s)])]  
group_al'' = sort . map (\(x,ys) -> (x,sort ys)) . groups . map (\(p,q,s) -> (p,(q,s))) . concat
-}



-- | Given a list of possible Alignments, for each p get entries of (q,s), such that q is unique. If not, take the q with the maximal s.
--   p,q,s values are sorted using an IntMap
--   As output, there is a vector of p's with a vector of (q,s) pairs for each p.
group_al''' :: V.Vector (Alignment s p q) -> V.Vector (Pos, U.Vector (Pos,Score))
group_al''' = list2vec . toList . M.unionsWith merge . map toMap . V.toList
    where merge = M.unionWith max -- (+) 
          toMap :: Alignment s p q -> M.IntMap (M.IntMap Score) -- =first coord (sec coor,sc)
          toMap = M.fromAscList . map (\(p,q,s) -> (p,M.singleton q s)) . (U.toList) . (U.map) (\(A p q s) -> (fromIntegral p,(fromIntegral q), s))
          toList :: M.IntMap (M.IntMap Score) -> [(Pos, [(Pos, Score)])]
          toList = map fst2int . M.toAscList . M.map (map fst2int . M.toAscList)
          fst2int (f,r) = (fromIntegral f,r)
          list2vec = V.fromList . map (\(a,b) -> (a, U.fromList b))





-- | the same as group_al''' but with a HashMap instead of an IntMap
group_al4 :: V.Vector (Alignment s p q) -> V.Vector (Pos, U.Vector (Pos,Score))
group_al4 = {-# SCC "group_al4_start" #-} list2vec . Data.List.sortBy (compare `on` fval) . toList . toIns . concatMap toMap . V.toList
    where toMap :: (Alignment s p q) -> [(Pos,(H.HashMap Pos Score))]
          toMap = {-# SCC "group_al4_toMap" #-} map (\(p,q,s) -> (p,H.singleton q s)) . (U.toList) . (U.map) (\(A p q s) -> (fromIntegral p,(fromIntegral q), s))
          toIns :: [(Pos,(H.HashMap Pos Score))] -> (H.HashMap Pos (H.HashMap Pos Score))
          toIns = {-# SCC "group_al4_toIns" #-} Data.List.foldl' (\a (p,v) -> H.insertWith f p v a) H.empty
          f = {-# SCC "group_al4_f" #-} H.unionWith max
          toList :: (H.HashMap Pos (H.HashMap Pos Score)) -> [(Pos,[(Pos,Score)])]
          toList = {-# SCC "group_al4_toList" #-} (H.toList . H.map (Data.List.sortBy (compare `on` fval) . H.toList))
          fst2int (f,r) = (fromIntegral f,r)
          list2vec = V.fromList . map (\(a,b) -> (a, U.fromList b))
          fval (f,_) = fromIntegral f
          f `on` g = \x y -> f (g x) (g y)
          


-- | the same as group_al''' but with a Judy array instead of an IntMap. p and q values are encoded as p*2^32 + q in the Judy array.
-- note:  (unsafePerformIO . withJ) or  unsafePerformIO . (\xs -> toJ xs >>= fromJ)
group_al5 :: V.Vector (Alignment s p q) -> V.Vector (Pos, U.Vector (Pos,Score))
group_al5 = {-# SCC "group_al5_start" #-} list2vec . editlist . unsafePerformIO . (\xs -> toJ xs >>= fromJ) . concat . map ( (U.toList) . (U.map) (\(A p q s) -> blubb (p,q,s)) ) . V.toList
    where
      list2vec :: [(Pos,[(Pos,Score)])] -> V.Vector (Pos, U.Vector (Pos,Score))
      list2vec = {-# SCC "group_al5_list2vec" #-} V.fromList . map (\(a,b) -> (a, U.fromList b)) 
      editlist :: [(Pos,[(Pos,Score)])] -> [(Pos,[(Pos,Score)])]
      editlist [] = []
      editlist [x] = [x]
      editlist zs = {-# SCC "group_al5_editlist" #-}  map (\g -> (fst $ head g, concatMap snd g)) $ groupBy  ((==) `on` fst) zs
      f `on` g = \x y -> f (g x) (g y)
      blubb :: (Pos,Pos,Score) -> (J.Key,Word32)
      blubb (p,q,s) = {-# SCC "group_al5_blubb" #-} ( (fromIntegral $ p `asTypeOf` p)*ck1 + (fromIntegral q), unsafeCoerce $ s )

ck1 :: Num a => a
ck1 = 2 ^ 32

-- | toJ and fromJ together in one function
withJ :: [(J.Key,Word32)] -> IO [(Pos, [(Pos,Score)])]
withJ xs = do
    -- performGC
    j <- J.new :: IO (J.JudyL Word32)
    mapM_ (\(k,s)-> do 
             v <- J.lookup k j
             case v of 
               Nothing -> J.insert k s j
               Just vv -> J.insert k (max s vv) j
          ) xs
    zs <- J.keys j >>= mapM (\k -> (format1 k . fromJust) <$> J.lookup k j) 
    return zs
        where
          format :: [J.Key] -> [Word32] -> [(Pos,[(Pos,Score)])] 
          format ks es = {-# SCC "fromJ_fromat" #-} map(\(k,v) -> let (p,q) = (convertback $ fromIntegral k)
                                                                  in (fromIntegral p,[(fromIntegral q, unsafeCoerce v ) ]) 
                                                       ) $ zip ks es

{-# INLINE convertback #-}
convertback :: Int -> (Int,Int)
convertback k = {-# SCC "fromJ_convertback" #-} quotRem k ck1

                                                                            
{-# INLINE format1 #-}                                                      
format1 :: J.Key -> Word32 -> (Pos,[(Pos,Score)])                           
format1 k v = let (p,q) = convertback $ fromIntegral k                      
              in  (fromIntegral p, [(fromIntegral q, unsafeCoerce v)])  

-- | convert a list into a Judy array
toJ xs = {-# SCC "toJ_start" #-} do
  performGC
  j <- J.new :: IO (J.JudyL Word32)
  mapM_ (\(k,s)-> do --k <- (pq)
           v <- J.lookup k j
           case v of 
             Nothing -> J.insert k s j
             Just vv -> J.insert k (max s vv) j
        ) xs
  return j

-- | make a list out of a Judy array
fromJ js = {-# SCC "fromJ_start" #-} do
  ks <- J.keys js
  -- es <- J.elems js 
  es <- catMaybes <$> sequence [ J.lookup k js | k <- ks ]
  return $ format ks es
      where
        format :: [J.Key] -> [Word32] -> [(Pos,[(Pos,Score)])] 
        format ks es = {-# SCC "fromJ_fromat" #-} map(\(k,v) -> let (p,q) = (convertback $ fromIntegral k)
                                                                   in (fromIntegral p,[(fromIntegral q, unsafeCoerce v ) ]) 
                                                        ) $ zip ks es
        convertback :: Int -> (Int,Int)
        convertback k = {-# SCC "fromJ_convertback" #-} quotRem k ck1


{-
a1, a2 :: Alignment Double Int Int
a1 = [(0,2,0.5),(1,3,0.4),(2,5,0.2)]
a2 = [(x,x+2,0.3) | x <- [0..3]]

test :: [(SeqId Int, Alignment Double Int Int)]
test = undefined "doit" [(undefined,a1)] (const [(undefined,a2),(undefined,a1)])

-- | Check colinearity and one-to-oneness by asserting that 
--   both sets of positions are sorted and without duplicates.
--   This won't necessarily work for alignments of a sequence to its reverse
--   complement, one way to do it is to represent a reverse match to position
--   k as (negate k-1), essentially laying out the sequence as:
--      <-----gcat-- ==atgc====>
--       -10 .... -1 0123456789
--   We still need to check that an alignment is either
--   all negative, or all non-negative.
prop_alignment :: (Ord p, Ord q) => Alignment s p q -> Bool
prop_alignment a = isSorted (map first a) && isSorted (map second a)
  where first (x,_,_) = x
        second (_,y,_) = y
        isSorted (x:y:rest) = x < y && isSorted (y:rest)
        isSorted _ = True
-}


--instance NFData B.ByteString where                                                                                                                         --  rnf a = a `seq` ()  

instance NFData A where
    rnf a = a `seq` ()
