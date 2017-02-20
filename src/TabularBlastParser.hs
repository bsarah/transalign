{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE LambdaCase #-}

-- This attoparsec module
module TabularBlastParser (module TabularBlastParser)
where
import Prelude hiding (takeWhile)
import Data.Attoparsec.ByteString.Char8 hiding (isSpace)
import qualified Data.ByteString.Char8 as C
import qualified Data.ByteString.Builder as S
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Vector as V
import System.Directory
import Data.Char

-- | reads and parses tabular Blast result from provided filePath
readTabularBlasts :: String -> IO (Either String (V.Vector BlastTabularResult))
readTabularBlasts filePath = do
  blastFileExists <- doesFileExist filePath
  if blastFileExists
     then do
       blastData <- C.readFile filePath
       let parsedBlast = parseTabularBlasts blastData
       return parsedBlast
     else return (Left ("Provided tabular blast file does not exist at:" ++ filePath))

parseTabularBlasts :: C.ByteString -> Either String (V.Vector BlastTabularResult)
parseTabularBlasts = parseOnly genParseTabularBlasts

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

genParseTabularBlasts :: Parser (V.Vector BlastTabularResult)
genParseTabularBlasts = do
  bresults <- many1 genParseTabularBlast
  return (V.fromList bresults)

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
  _tabularHit <- many' (try genParseBlastTabularHit)  <?> "Tabular hit"
  return $ BlastTabularResult _blastProgram (toLB $ _blastQueryId) (toLB $ C.pack _blastDatabase) _blastHitNumber (V.fromList _tabularHit)

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
  _queryId <- many1 (notChar '\t') <?> "hit qid"
  char '\t'
  _subjectId <- many1 (notChar '\t') <?> "hit sid"
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
  _querySeq <- many1 (satisfy bioLetters) <?> "hit qseq"
  char '\t'
  _subjectSeq <- many1 (satisfy bioLetters) <?> "hit subSeq"
  char '\n'
  return $ BlastTabularHit (B.pack _queryId) (B.pack _subjectId) _seqIdentity _alignmentLength _misMatches _gapOpenScore _queryStart _queryEnd _hitSeqStart _hitSeqEnd _eValue _bitScore _subjectFrame (B.pack _querySeq) (B.pack _subjectSeq)
  
--IUPAC amino acid with gap
aminoacidLetters :: Char -> Bool
aminoacidLetters = inClass "ARNDCQEGHILMFPSTWYVBZX-"

--IUPAC nucleic acid characters with gap
nucleotideLetters :: Char -> Bool
nucleotideLetters = inClass "AGTCURYSWKMBDHVN-."

--IUPAC nucleic acid characters with gap
bioLetters :: Char -> Bool
bioLetters = inClass "ABCDEFGHIJKLMNOPQRSTVWXYZ.-"


toLB :: C.ByteString -> B.ByteString
toLB = S.toLazyByteString . S.byteString
