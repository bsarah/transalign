{-# Language OverloadedStrings #-}
{-# Language TypeFamilies #-}
{-# Language MultiParamTypeClasses #-}
{-# Language TemplateHaskell #-}
{-# Language FlexibleInstances #-}
-- {-# Language UndecidableInstances #-}                                                                                                                     
-- {-# Language FlexibleContexts #-}                                                                                                                         
{-# LANGUAGE CPP #-}


module Blast where

import qualified Data.ByteString.Lazy.Char8 as B
import Bio.BlastXML -- hiding (SeqId)
import Bio.Core
import qualified Data.Map as M
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import Data.Vector.Unboxed.Deriving
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable


import Data.Int
import Text.Printf
import GHC.Float (double2Float)
import Align


type BlastAlignData = Alignment Float Int32 Int32
type BlastAlignment = (B.ByteString,BlastAlignData)

newtype SeqId a = S a deriving (Show,Ord,Eq)

qname :: BlastRecord -> B.ByteString
qname = head . B.words . unSL . query

sname :: BlastHit -> B.ByteString
sname = head . B.words . unSL . subject

data BlastProg = BlastX | BlastP -- etc

hit2align :: BlastProg -> BlastHit -> [BlastAlignment]
hit2align p h = [(sname h, a) | a <- map (match2align p) $ matches h]
        
-- | Construct an alignment from a BlastMatch data structure.
-- Warning: blastx-specific!
-- BlastMatch has two attributes qseq and hseq :: ByteString, additionally q_from,q_to,h_from,h_to as Int 
match2align :: BlastProg -> BlastMatch -> BlastAlignData
match2align prog m = G.map add_score $ go qstart hstart (qseq m) (hseq m)
  where (hstep,qstep) = case prog of BlastX -> (1,3) -- blastx? (1,1) else
                                     BlastP -> (1,1)
                                     -- etc.
        go qi hi qsq hsq = case (B.uncons qsq,B.uncons hsq) of
          (Just ('-',qs),Just (_h,hs)) -> go qi (hi+hstep) qs hs
          (Just (_q,qs),Just ('-',hs)) -> go (qi+qstep) hi qs hs
          (Just (_q,qs),Just (_h,hs))   -> (qi,hi) `G.cons` go (qi+qstep) (hi+hstep) qs hs
          (Nothing,Nothing)           -> G.empty
          _ -> error "Case closed! just to shut up the type checker"
        qstart = case aux m of Frame Minus _ -> negate (q_to m)
                               _ -> q_from m
        hstart = h_from m
        hit_length = (h_to m- h_from m + 1)
        add_score (p,q) = (A (fromIntegral p) (fromIntegral q) (double2Float (bits m/fromIntegral hit_length)))
        
{-
-- merge different HSPs
-- each position in query must map to a unique position in target        
-- in case of ambiguity, pick best hit, but weaken it
-- DO NOT USE!        
unify :: (Ord q, Ord p, Ord s, Num s) => [Alignment s p q] -> Alignment s p q
unify = uniq . sort . concat
  where uniq ((x,y,z):(a,b,c):rest) 
          | x==a = if z < c then uniq ((a,b,c-z):rest)
                   else uniq ((x,y,z-c):rest)
          | otherwise = (x,y,z):uniq ((a,b,c):rest)
        uniq [x] = [x]
        uniq [] = []
-}
          
--                     seqid             seqid       alignment
type BlastMap = M.Map B.ByteString [BlastAlignment]
type BlastId = B.ByteString

targets :: BlastMap -> [SeqId B.ByteString]
targets = map (S . fst) . concat . map snd . M.toList

-- | Extract the set of alignments from a Blast XML file
getAlignments :: BlastResult -> [(B.ByteString, [BlastAlignment])]
getAlignments res = map rec2align . results $ res
  where rec2align r = (qname r,concatMap (hit2align prog) $ hits r)
        prog = case blastprogram res of "blastp" -> BlastP
                                        "blastx" -> BlastX
                                        _ -> error ("undefined blastprogram")

-- | Read a set of alignments from a BlastXML file
readAlignments :: FilePath -> IO [(B.ByteString, [BlastAlignment])]
readAlignments f = getAlignments `fmap` readXML f

-- | Read alignments, and return a Map for query to set of alignments
readAlignmentMap :: FilePath -> IO BlastMap
readAlignmentMap f = M.fromList `fmap` getAlignments `fmap` readXML f

-- | Read Blast alignments, but only for the given seqids, and with a max no of hits per seqid
readFilteredAlignmentMap :: [SeqId B.ByteString] -> Int -> FilePath -> IO BlastMap
readFilteredAlignmentMap seqs maxnum f = do
  M.fromList `fmap` extractBlast `fmap` readXML f
  where extractBlast res = map (rec2align (prog res)) . filter wanted . results $ res
        wanted r = S (qname r) `elem` seqs
        rec2align p r = (qname r,concatMap (hit2align p) $ take maxnum $ hits r)
        prog r = case blastprogram r of "blastp" -> BlastP
                                        "blastx" -> BlastX
                                        _ -> error ("undefined blastprogram")

-- testing:

{-
test0 :: IO [BlastRecord]
test0 = results `fmap` readXML "contigs_vs_sprot.xml"

test1 :: BlastRecord -> Maybe BlastRecord
test1 xs = if any not . map (prop_alignment . snd) . concatMap (hit2align BlastX) . hits $ xs
           then Just xs else Nothing

test2 :: IO BlastMap
test2 = readAlignmentMap "contigs_vs_sprot.xml"
-}

print_alignment :: B.ByteString -> BlastAlignment -> IO ()
print_alignment src (tgt,al) = do
  printf "%s\t%s\n" (B.unpack src) (B.unpack tgt)
  G.mapM_ (printf "%6d\t" . fst') al
  printf "\n"
  G.mapM_ (printf "%6.2f\t" . trd') al
  printf "\n"
  G.mapM_ (printf "%6d\t" . snd') al
  printf "\n"
  
fst' :: A -> Int32
fst' (A x _ _) = x

snd' :: A -> Int32
snd' (A _ y _) = y

trd' :: A -> Float
trd' (A _ _ z) = z
