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
{-# LANGUAGE CPP #-}

module Align where

import qualified Data.IntMap.Strict as M
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import Data.Vector.Unboxed.Deriving
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable
import Data.List (sort,sortBy)
import Data.Int (Int32)
import Debug.Trace
import Control.DeepSeq
import Data.Vector.Algorithms.Merge


type Pos = Int32
type Score = Float

data A  = A {-# UNPACK #-} !Pos  {-# UNPACK #-} !Pos  {-# UNPACK #-} !Score deriving Show

derivingUnbox "A"
    [t| A -> (Pos,Pos,Score) |]
    [| \(A p q s) -> (p,q,s) |]
    [| \(p,q,s) -> A p q s |]

type Alignment score pos1 pos2 = U.Vector A


collect_aligns :: (usid -> [(ssid, Alignment s q r)]) -> [(usid, Alignment s p q)] -> [(ssid,Alignment s p r)]
collect_aligns sp_lu ups = do -- list monad
  (uhit,ualign) <- ups
  (shit,salign) <- sp_lu uhit
  return (shit,trans_align ualign salign)

-- A proper 'Alignment' links positions in one sequence to positions in another,
-- assigning a score to each link.  The position links are co-linear and one-to-one.


score :: Alignment s p q -> Score -- s
score = G.sum . G.map (\(A _ _ s) -> s)

-- | Reverse the alignment
invert :: Alignment s p q -> Alignment s q p
invert = G.map (\(A x y s) -> (A y x s))

-- | Follow positions via the first 'Alignment', through the second to 
--   form the transitive alignment.  Note that if the alignments have 
--   nothing in common, this may produce an empty alignment.
{-
trans_align :: Alignment s p q -> Alignment s q r -> Alignment s p r
trans_align x@((A xp xq xs):xx) y@((A yq yr ys):yy)
  | xq == yq = (A xp yr (min xs ys)) : trans_align xx yy
  | xq <  yq = trans_align xx y
  | xq >  yq = trans_align x yy
trans_align _ _ = [] -}
-- NOTE: this is exactly the >>> operator for the "Alignment s" arrow!

{-
trans_align :: Alignment s p q -> Alignment s q r -> Alignment s p r
trans_align x y 
  | G.null x || G.null y = G.empty
  | otherwise = let (A xp xq xs) = G.head x
                    (A yq yr ys) = G.head y
                    xt = G.tail x
                    yt = G.tail y
                in
                 if (xq == yq) then (A xp yr (min xs ys))  `G.cons` trans_align xt yt
                 else 
                   if (xq < yq) then trans_align xt y
                   else {-xq > yq =-} trans_align x yt

-}


trans_align :: Alignment s p q -> Alignment s q r -> Alignment s p r
trans_align xorig yorig = G.unfoldr go (xorig,yorig) where
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





-- collect hits against the same sequence, and calculate a consensus alignment
merge_aligns :: (Show ssid, Ord ssid) => [(ssid,Alignment s p q)] -> [(ssid,Alignment s p q)]
merge_aligns = map merge_group . groups . filter is_not_empty --groups:: [(sid,align)] -> [(sid,V.Vector align)]
  where is_not_empty = not . G.null . snd



merge_group :: (Show ssid) => (ssid, V.Vector (Alignment s p q)) -> (ssid, Alignment s p q)
merge_group (tgt,as) = let bs = group_al''' as 
                       in (tgt , go (V.fromList []) bs) -- $ map (\(p,xs) -> (p,collect xs)) 
    where        
    -- this is just a regular alignment using the 
    -- scores from the provided alignments. No sparse matrix.
   -- go :: V.Vector (score, Alignment s p q) -> V.Vector (Pos, V.Vector (Pos,Score)) -> Alignment s p q
    go y x
      | G.null y && G.null x = error ("last will fail! " ++ show tgt ++ " " ++ show as)
      | G.null x = U.reverse $ snd $ V.last $ y
      | G.null y = let (p,xs) = G.head x
                       rest = G.tail x
                   in go (G.map (\(q,s)->(s,G.singleton (A p q s))) xs) rest                           
      | otherwise = let xh = G.head x    
                        xt = G.tail x
                    in go (sort_column $ merge1 y xh) xt
        
    

-- Merge two columns (the previous column and the next one) in the alignment DP matrix                                                                      
-- (t ~ Pos, p ~ Pos, q ~ Pos, a ~ Pos, s ~ Score, t1 ~ Score,  Num t1, Ord t1, Ord a) =>                                                                   
-- warning: invariant: total-score ts increases down the prev_column!                                                                                        
merge1 :: V.Vector (Score, U.Vector A) -> (Pos, V.Vector (Pos, Score)) -> V.Vector (Score, U.Vector A) --dynamic programming, forward tracking, not many cos\ts:      
merge1 x y = V.unfoldr go (Just(x,y))
    where go Nothing = Nothing --go :: Maybe(x,y) -> Maybe((score,uvec A),Maybe(x,y))                                                                       
          go (Just(x,y))
              | V.null $ snd y = if V.null x then Nothing else Just(V.head x,Just(V.tail x,y))
              | V.null x = error "error in merge1, empty prev"
              | G.null $ snd $ V.head x = error "merge1: complex pattern missing"
              | q <= q1 = Just((s,U.singleton (A p q s)),Just( x,(p,nc))) --prev@((ts,(A p1 q1 s1):qs1):pc) (p,(q,s):nc), x=prev                            
              | V.null pc = Just( (ts + s, G.unfoldr go' (Just((A p q s),as)) ) , Just( x, (p,nc)) )
              | otherwise = let (_,hpc) = V.head pc
                                (A _ q2 _) = U.head hpc
                            in if (q2 < q)
                               then go (Just(pc,(p,ys)))
                               else Just((ts + s,(G.unfoldr go' (Just((A p q s),as)))) ,Just( x, (p,nc)))
             where (ts,as) = V.unsafeHead x
                    pc = V.unsafeTail x
                    (A p1 q1 s1) = U.head as
                    qs1 = U.tail as
                    (p,ys) = y
                    (q,s) = V.head ys
                    nc = V.tail ys
                    go' Nothing = Nothing
                    go' (Just(aa,bb))
                        |G.null bb = Just(aa,Nothing)
                        |otherwise = Just(aa,Just(G.head bb,G.tail bb))




-- sort on q-coordinate and filter so scores are increasing                                                                                                 
sort_column :: V.Vector (Score, U.Vector A) -> V.Vector (Score, U.Vector A)
sort_column = filter_scores . sortBy' (compare `on` qval)
  where qval (sc,as) = let (A p q s) = U.head as
                       in q
        f `on` g = \x y -> f (g x) (g y)
        sortBy' f xs = G.modify (Data.Vector.Algorithms.Merge.sortBy f) xs 
        filter_scores x = if V.null x then V.empty else V.unfoldr go (Just(V.head x,V.tail x))
            where go Nothing = Nothing
                  go (Just(a,b))
                      | V.null b = Just(a,Nothing)
                      | otherwise = let (xs,xss) = a
                                        (ys,yss) = V.head b
                                        zzs = V.tail b
                                    in if ys < xs then go (Just((xs,xss),zzs))
                                       else Just((xs,xss),Just((ys,yss),zzs))
    

        
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
 -- merge :: Alignment s p q -> Alignment s p q -> Alignment s p q
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


-}

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
group_al''' :: V.Vector (Alignment s p q) -> V.Vector (Pos, V.Vector (Pos,Score))
--group_al''' :: [Alignment s p q] -> [(Pos,[(Pos,Score)])] --dp matrix, first pos = column
--group_al''' xs = let y = list2vec . toList . M.unionsWith merge . map toMap . V.toList $ xs in traceShow("groupal", V.length y) y
group_al'''= list2vec . toList . M.unionsWith merge . map toMap . V.toList
    where merge = M.unionWith max -- (+) 
          toMap :: Alignment s p q -> M.IntMap (M.IntMap Score) -- =first coord (sec coor,sc)
          toMap = M.fromAscList . map (\(p,q,s) -> (p,M.singleton q s)) . (U.toList) . (U.map) (\(A p q s) -> (fromIntegral p,(fromIntegral q), s))
          toList :: M.IntMap (M.IntMap Score) -> [(Pos, [(Pos, Score)])]
          toList = map fst2int . M.toAscList . M.map (map fst2int . M.toAscList)
          fst2int (f,r) = (fromIntegral f,r)
          list2vec = V.fromList . map (\(a,b) -> (a, V.fromList b))



{-
group_al''' :: V.Vector (Alignment s p q) -> V.Vector (Pos, V.Vector (Pos,Score))
group_al''' xs = let y = foldi . sortBy' (compare `on` pval) .  G.concat . G.toList . V.map toTup . V.map  toTup' $ xs in traceShow("groupal", V.length y) y
    where pval (p,_) = p
          toTup' :: (Alignment s p q) -> V.Vector(Pos,Pos,Score)
          toTup' =  V.fromList . U.toList . U.map (\(A p q s) -> ((fromIntegral p), (fromIntegral q), s))
          toTup :: V.Vector(Pos,Pos,Score) -> V.Vector(Pos,V.Vector((Pos,Score)))
          toTup = V.map (\(p,q,s) -> (p, V.singleton(q,s)))
          f `on` g = \x y -> f (g x) (g y)
          sortBy' f xs = G.modify (Data.Vector.Algorithms.Merge.sortBy f) xs --(V.fromList) $ sortBy f (V.toList xs)                                        
          foldi xs = V.unfoldr go xs
              where
                go xs
                    | V.length xs == 1 = let (p1,vs1) = V.head xs
                                             txs = V.tail xs
                                         in Just((p1,vs1),txs2)
                    | V.null xs = Nothing
                    | p1 == p2 = Just((p1,(V.unfoldr go'(vs1,vs2))), t2)
                    | otherwise = Just((p1,vs1),txs2)
                    where (p1,vs1) = V.head xs
                          txs2 = V.tail xs
                          (p2,vs2) = V.head txs2
                          t2 = V.tail txs2
                          go' (vs1,vs2)
                              | V.null vs1 && V.null vs2 = Nothing
                              | V.null vs1 = Just(V.head vs2,(V.tail vs2,vs1))
                              | otherwise = Just(V.head vs1,(V.tail vs1,vs2))
-}







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
