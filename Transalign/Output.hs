module Transalign.Output where

import qualified Data.ByteString.Lazy.Char8 as B
import Text.Printf (printf)
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import Data.Vector.Unboxed.Deriving
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable


import Transalign.Blast (BlastAlignData, print_alignment)
import Transalign.Align (Alignment,score,A(..))

output_long :: B.ByteString -> [(Float,B.ByteString,BlastAlignData)] -> IO ()
output_long src = do 
  mapM_ (\(s,sid,xs) -> print_alignment src (sid,xs) >> printf "score: %f\n" s)

output_short :: B.ByteString -> [(Float,B.ByteString,BlastAlignData)] -> IO ()
output_short src = mapM_ (print1 src)

print1 :: B.ByteString -> (Float, B.ByteString, BlastAlignData) -> IO ()
print1 src (sc,tgt,al) = do
  let (A p0 q0 _s0) = G.head al
      (A pn qn _sn) = G.last al
  printf "%s\t%s\t%8.2f\t %6d\t %3.2f\t %5d\t %5d\t %5d\t %5d\n"
    (B.unpack src) (B.unpack tgt) sc (G.length al) (sc/fromIntegral (G.length al)) p0 pn q0 qn
  -- todo: global scores (div by protein length and transcript length)
  -- and local score (div by alignment length)
  
-- score alignments
add_score :: (sid,Alignment t p q) -> (Float, sid, Alignment t p q)
add_score (s,x) = (score x, s, x)
-- prop: score . hit2align == score blasthit
