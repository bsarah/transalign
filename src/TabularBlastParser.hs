{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}

-- This attoparsec module
module TabularBlastParser (module TabularBlastParser)
where
import Prelude hiding (takeWhile)
import Data.Attoparsec.ByteString.Char8 hiding (isSpace)
import qualified Data.Attoparsec.ByteString.Lazy as L
import qualified Data.ByteString.Char8 as C
import qualified Data.ByteString.Builder as S
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Vector as V
import System.Directory
import Data.Char
import Control.Monad
import Debug.Trace
import Text.Printf

-- | reads and parses tabular Blast result from provided filePath
readTabularBlasts :: String -> IO [BlastTabularResult]
readTabularBlasts filePath = do
  printf "# reading tabular blast input from file %s\n" filePath
  blastFileExists <- doesFileExist filePath
  if blastFileExists
     then parseTabularBlasts <$> B.readFile filePath
     else fail "# tabular blast file \"%s\" does not exist\n" filePath

-- |Read a lazy bytestring and stream out a lsit of @BlastTabularResult@'s.
-- In case, there is a parse error "late" in the file, we might have
-- already streamed out some (or many!) of these results.

parseTabularBlasts :: B.ByteString -> [BlastTabularResult]
parseTabularBlasts = go
  where go xs = case L.parse genParseTabularBlast xs of
          L.Fail remainingInput ctxts err  -> error $ "parseTabularBlasts failed! " ++ err ++ " ctxt: " ++ show ctxts ++ " head of remaining input: " ++ (B.unpack $ B.take 1000 remainingInput)
          L.Done remainingInput btr ->
            btr : go remainingInput

data BlastTabularResult = BlastTabularResult
  { blastProgram :: !BlastProgram,
    blastQueryId :: !B.ByteString,
--    blastQueryName :: !B.ByteString,
    blastDatabase :: !B.ByteString,
    blastHitNumber :: !Int,
    hitLines :: !(V.Vector BlastTabularHit)
  }
  deriving (Show, Eq)

data BlastProgram = BlastX | BlastP | BlastN
  deriving (Show, Eq)

genParseTabularBlasts :: Parser [BlastTabularResult]
genParseTabularBlasts = do
  bresults <- many1 genParseTabularBlast
  return bresults

genParseBlastProgram :: Parser BlastProgram
genParseBlastProgram = do
  choice [string "# BLAST",string "# blast"]
  (toLower <$> anyChar) >>= return . \case 
    'x' -> BlastX
    'p' -> BlastP
    'n' -> BlastN

genParseTabularBlast :: Parser BlastTabularResult
genParseTabularBlast = do
  --choice [string "# BLAST",string "# blast"]
  --_blastProgram <- choice [string "p",string "P",string "x",string "X"] <?> "Program"
  _blastProgram <- genParseBlastProgram <?> "Program"
  many1 (notChar '\n')
  endOfLine
  string "# Query: " <?> "Query"
  --_blastQueryId <- many1 (notChar ' ') <* skipWhile (not (string "\n")) <?> "QueryId"
  _blastQueryId <- takeWhile (not . isSpace) <* manyTill anyChar endOfLine <?> "QueryId"
  --_blastQueryName <- many' (notChar '\n')  <?> "QueryName"
  string "# Database: " <?> "Database"
  _blastDatabase <- many1 (notChar '\n') <?> "Db"
  string "\n# " <?> "header linebreak"
  --fields line
  skipMany (try genParseFieldLine) <?> "Fields"
  _blastHitNumber <- decimal  <?> "Hit number"
  string " hits found\n" <?> "hits found"
  _tabularHit <- count  _blastHitNumber (try genParseBlastTabularHit)  <?> "Tabular hit"
  return $ BlastTabularResult _blastProgram (toLB $ _blastQueryId) (B.pack _blastDatabase) _blastHitNumber (V.fromList _tabularHit)

genParseFieldLine :: Parser ()
genParseFieldLine = do
  string "Fields:"
  skipMany (notChar '\n')
  string "\n# "
  return ()

data BlastTabularHit = BlastTabularHit
  { queryId :: !B.ByteString,
    subjectId ::  !B.ByteString,
    seqIdentity :: !Double,
    alignmentLength :: !Int,
    misMatches :: !Int,
    gapOpenScore :: !Int,
    queryStart :: !Int,
    queryEnd :: !Int,
    hitSeqStart :: !Int,
    hitSeqEnd :: !Int,
    eValue :: !Double,
    bitScore :: !Double,
    subjectFrame :: !Int,
    querySeq  :: !B.ByteString,
    subjectSeq  :: !B.ByteString
  }
  deriving (Show, Eq)

genParseBlastTabularHit :: Parser BlastTabularHit
genParseBlastTabularHit = do
  --start <- peekChar'
  --when (start == '#') (fail "irgendwas")
  --_queryId <- takeWhile (\c ->  c /= '#' || c /= '\t')  <?> "hit qid"
  _queryId <- takeWhile1 ((/=9) . ord) <?> "hit qid"
  char '\t'
  _subjectId <- takeWhile1 ((/=9) . ord) <?> "hit sid"
  char '\t'
  _seqIdentity <- double <?> "hit seqid"
  char '\t'
  _alignmentLength <- decimal  <?> "hit sid"
  char '\t'
  _misMatches <- decimal <?> "hit mmatch"
  char '\t'
  _gapOpenScore <- decimal <?> "hit gopen"
  char '\t'
  _queryStart <- decimal <?> "hit qstart"
  char '\t'
  _queryEnd <- decimal  <?> "hit qend"
  char '\t'
  _hitSeqStart <- decimal  <?> "hit sstart"
  char '\t'
  _hitSeqEnd <- decimal <?> "hit send"
  char '\t'
  _eValue <- double <?> "hit eval"
  char '\t'
  _bitScore <- double <?> "hit bs"
  char '\t'
  _subjectFrame <- decimal <?> "hit sF"
  char '\t'
  _querySeq <- takeWhile1 ((/=9) . ord) <?> "hit qseq" -- 9 == '\t'
  char '\t'
  _subjectSeq <- takeWhile1 ((/=10) . ord) <?> "hit subSeq" -- 10 == '\n'
  char '\n'
  return $ BlastTabularHit (B.fromStrict _queryId) (B.fromStrict _subjectId) _seqIdentity _alignmentLength _misMatches _gapOpenScore _queryStart _queryEnd _hitSeqStart _hitSeqEnd _eValue _bitScore _subjectFrame (B.fromStrict _querySeq) (B.fromStrict _subjectSeq)
  
--IUPAC amino acid with gap
--aminoacidLetters :: Char -> Bool
aminoacidLetters = inClass "ARNDCQEGHILMFPSTWYVBZX-"

--IUPAC nucleic acid characters with gap
--nucleotideLetters :: Char -> Bool
nucleotideLetters = inClass "AGTCURYSWKMBDHVN-."

--IUPAC nucleic acid characters with gap
--bioLetters :: Char -> Bool
bioLetters = inClass "ABCDEFGHIJKLMNOPQRSTUVWXYZ.-"


toLB :: C.ByteString -> B.ByteString
toLB = S.toLazyByteString . S.byteString
