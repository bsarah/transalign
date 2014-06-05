{- IDEA:  annotate a transcript by doing blastx vs NR, then blastp NR vs SwissProt
   Align transcript vs SP by mapping each position via NR hits to positions in SP.
   How to test?

   Generalize: 
     annotate a_vs_b.xml b_vs_c.xml c_vs_d.xml => result: d_vs_a (or: a_vs_d?)
-}

{-# Language TypeFamilies #-}

module Align where

import qualified Data.IntMap.Strict as M
import Data.List (sort,sortBy)
import Data.Int (Int32)

type Pos = Int32
type Score = Float

collect_aligns :: (usid -> [(ssid, Alignment s q r)]) -> [(usid, Alignment s p q)] -> [(ssid,Alignment s p r)]
collect_aligns sp_lu ups = do -- list monad
  (uhit,ualign) <- ups
  (shit,salign) <- sp_lu uhit
  return (shit,trans_align ualign salign)

-- A proper 'Alignment' links positions in one sequence to positions in another,
-- assigning a score to each link.  The position links are co-linear and one-to-one.

data A = A {-# UNPACK #-} !Pos  {-# UNPACK #-} !Pos  {-# UNPACK #-} !Score deriving Show
type Alignment score pos1 pos2 = [A]

score :: Alignment s p q -> Score -- s
score = sum . map (\(A _ _ s) -> s)

-- | Reverse the alignment
invert :: Alignment s p q -> Alignment s q p
invert = map (\(A x y s) -> (A y x s))

-- | Follow positions via the first 'Alignment', through the second to 
--   form the transitive alignment.  Note that if the alignments have 
--   nothing in common, this may produce an empty alignment.
trans_align :: Alignment s p q -> Alignment s q r -> Alignment s p r
trans_align x@((A xp xq xs):xx) y@((A yq yr ys):yy)
  | xq == yq = (A xp yr (min xs ys)) : trans_align xx yy
  | xq <  yq = trans_align xx y
  | xq >  yq = trans_align x yy
trans_align _ _ = []
-- NOTE: this is exactly the >>> operator for the "Alignment s" arrow!

-- collect hits against the same sequence, and calculate a consensus alignment
merge_aligns :: (Show ssid, Ord ssid) => [(ssid,Alignment s p q)] -> [(ssid,Alignment s p q)]
merge_aligns = map merge_group . groups . filter is_not_empty
  where is_not_empty = not . null . snd

merge_group :: (Show ssid) => (ssid, [Alignment s p q]) -> (ssid, Alignment s p q)
merge_group (tgt,as) = (tgt , go [] $ group_al''' as) -- $ map (\(p,xs) -> (p,collect xs)) 
  where        
    -- this is just a regular local alignment using the 
    -- scores from the provided blast alignments.  Using the (sparse) matrix from group_al,
    -- no gap penalties, just the best increasing path
    go [] [] = error ("last will fail! " ++ show tgt ++ " " ++ show as) 
    go [] ((p,xs):rest) = go (map (\(q,s)->(s,[(A p q s)])) xs) rest     -- build the first columnt
    go prev_col (qss:rest) = go (sort_column $ merge1 prev_col qss) rest -- merge, and sort by q - i.e. target position
    go prev_col  [] = reverse $ snd $ last $ prev_col                    -- done: return the resulting alignment
    
-- Merge two columns (the accumulated column and the next one) from the alignment DP matrix
-- The accumulated column is a list of (scores,partial alignment)-pairs.  Each pair (cell) from the next
-- column is added to the highest scoring compatible entry in the accumulating column.
-- (t ~ Pos, p ~ Pos, q ~ Pos, a ~ Pos, s ~ Score, t1 ~ Score,  Num t1, Ord t1, Ord a) => 
merge1 :: [(Score, [A])] -> (Pos, [(Pos, Score)]) -> [(Score, [A])]
merge1 [] _ = error "empty prev in merge1"
merge1 ((_, []) : _) (_, _ : _) = error "merge1: complex pattern missing"
merge1 ps (_,[]) = ps 
merge1 prev@((ts,(A p1 q1 s1):qs1):pc) (p,(q,s):nc)
      | q <= q1 = (s,[(A p q s)]) : continue         -- q is the shortest alignment, and we just add this
      | null pc = (ts+s,(A p q s):(A p1 q1 s1):qs1) : continue  -- we exhausted prev, tack this pair on to the last entry
      | otherwise = let (_,(A _ q2 _):_) = head pc  -- now: let's see
                    in if q2 < q {- && ts2 > ts -} 
                       then merge1 pc (p,(q,s):nc)  -- prune this position, we're better off with the next one(?)
                       else (ts+s,(A p q s):(A p1 q1 s1):qs1) : continue  -- this is the best candidate, add it here
  where continue = merge1 prev (p,nc) -- recurse (so this is really just a fold)

-- warning: required invariant: total-score ts increases down the prev_column! (because the 
-- q-position increases, and for any alignment to be useful, it must score higher than shorter 
-- alternatives.)

-- sort on q-coordinate and filter so scores are increasing
sort_column :: [(Score, [A])] -> [(Score, [A])]
sort_column = filter_scores . sortBy (compare `on` qval)
  where qval (_ts,(A _p q _s):_qs) = q
        f `on` g = \x y -> f (g x) (g y)
        filter_scores ((xs,xss):(ys,yss):rest)
          | ys < xs = filter_scores ((xs,xss):rest)
          | otherwise = (xs,xss):filter_scores ((ys,yss):rest)
        filter_scores [x] = [x]
        filter_scores []  = []
          
-- Add scores for all (q,s) pairs mapping to the same q. Input must be sorted on q.
collect :: [(Pos, Score)] -> [(Pos, Score)]
collect ((q1,s1):(q2,s2):rest)  
  | q1==q2    = collect ((q1,s1+s2):rest)
  | q1 < q2   = (q1,s1) : collect ((q2,s2):rest)
  | otherwise = error ("collect: unsorted input!\n"++show ((q1,s1):(q2,s2):take 2 rest))
collect [x] = [x]
collect []  = []


-- Group alignments against the same target
-- Basically a quicksort
groups :: Ord sid => [(sid, align)] -> [(sid,[align])]
groups ((s,a):xs) = this : groups less ++ groups more
    where this = (s, a : map snd (filter ((s==).fst) xs))
          less = filter ((s<).fst) xs
          more = filter ((s>).fst) xs
groups [] = []
-- prop: length = length . concat . groups

-- Group alignments on first coordinate (p), using a direct merge of
-- the alignments.  In other words, we take a list of input alignments, and populate a
-- a complete, sparse, DP matrix.  (The group_alX functions, with X being a number 
-- of primes, are semantically equivalent implementations.
group_al :: [Alignment s p q] -> [(Pos,[(Pos,Score)])]
group_al = join . foldr1 merge
  where -- merge :: Alignment s p q -> Alignment s p q -> Alignment s p q
        merge one@((A p1 q1 s1):one_rest) two@((A p2 q2 s2):two_rest) 
          | p1 < p2   = (A p1 q1 s1):merge one_rest two
          | p1 == p2  = (A p1 q1 s1):(A p2 q2 s2):merge one_rest two_rest
          | p1 > p2 = (A p2 q2 s2):merge one two_rest
        merge [] xs = xs
        merge xs [] = xs
        join ((A p q s):xs) = go p [(q,s)] xs
        go p0 acc ((A p q s):xs) 
          | p0 == p   = go p0 ((q,s):acc) xs
          | otherwise = (p0,sort acc) : join ((A p q s):xs)
        go p0 acc [] = [(p0,sort acc)]
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

-- Use a Map p (Map q s): this is faster and uses less memory than the others.
group_al''' :: [Alignment s p q] -> [(Pos,[(Pos,Score)])]
group_al''' = toList . M.unionsWith merge . map toMap
  where merge = M.unionWith max -- (+)
        toMap :: Alignment s p q -> M.IntMap (M.IntMap Score)
        toMap = M.fromAscList . map (\(A p q s) -> (fromIntegral p,M.singleton (fromIntegral q) s))
        toList :: M.IntMap (M.IntMap Score) -> [(Pos, [(Pos, Score)])]
        toList = map fst2int . M.toAscList . M.map (map fst2int . M.toAscList)
        fst2int (f,r) = (fromIntegral f,r)
-- testing

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
