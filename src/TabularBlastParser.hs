{-# LANGUAGE OverloadedStrings #-}

-- This attoparsec module
module TabularBlastParser (module TabularBlastParser)
where
import Data.Attoparsec.ByteString.Char8
import qualified Data.ByteString.Char8 as C
import qualified Data.Vector as V
import System.Directory

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
  { blastProgram :: C.ByteString,
    blastQuery :: C.ByteString,
    blastDatabase :: C.ByteString,
    blastHitNumber :: Int,
    hitLines :: V.Vector BlastTabularHit
  }
  deriving (Show, Eq)

genParseTabularBlasts :: Parser (V.Vector BlastTabularResult)
genParseTabularBlasts = do
  bresults <- many1 genParseTabularBlast
  return (V.fromList bresults)

genParseTabularBlast :: Parser BlastTabularResult
genParseTabularBlast = do
  choice [string "# BLAST",string "# blast"]
  _blastProgram <- many1 (notChar '\n') <?> "Program"
  string "\n# Query: "
  _blastQuery <- many1 (notChar '\n')  <?> "Query"
  string "\n# Database: "
  _blastDatabase <- many1 (notChar '\n')  <?> "Db"
  string "\n# "
  --fields line
  skipMany (try genParseFieldLine)
  _blastHitNumber <- decimal  <?> "Hit number"
  string " hits found\n"
  _tabularHit <- many' genParseBlastTabularHit
  return $ BlastTabularResult (C.pack ("BLAST" ++ _blastProgram)) (C.pack _blastQuery) (C.pack _blastDatabase) _blastHitNumber (V.fromList _tabularHit)

genParseFieldLine :: Parser ()
genParseFieldLine = do
  string "Fields:"
  skipMany (notChar '\n')
  string "\n# "
  return ()

data BlastTabularHit = BlastTabularHit
  { queryId :: C.ByteString,
    subjectId ::  C.ByteString,
    seqIdentity :: Double,
    alignmentLength :: Int,
    misMatches :: Int,
    gapOpenScore :: Int,
    queryStart :: Int,
    queryEnd :: Int,
    hitSeqStart :: Int,
    hitSeqEnd :: Int,
    eValue :: Double,
    bitScore :: Double,
    subjectFrame :: Int,
    querySeq  :: C.ByteString,
    subjectSeq  :: C.ByteString
  }
  deriving (Show, Eq)

genParseBlastTabularHit :: Parser BlastTabularHit
genParseBlastTabularHit = do
  _queryId <- many1 (notChar '\t')
  char '\t'
  _subjectId <- many1 (notChar '\t')
  char '\t'
  _seqIdentity <- double
  char '\t'
  _alignmentLength <- decimal
  char '\t'
  _misMatches <- decimal
  char '\t'
  _gapOpenScore <- decimal
  char '\t'
  _queryStart <- decimal
  char '\t'
  _queryEnd <- decimal
  char '\t'
  _hitSeqStart <- decimal
  char '\t'
  _hitSeqEnd <- decimal
  char '\t'
  _eValue <- double
  char '\t'
  _bitScore <- double
  char '\t'
  _subjectFrame <- decimal
  char '\t'
  _querySeq <- many1 (satisfy bioLetters)
  char '\t'
  _subjectSeq <- many1 (satisfy bioLetters)
  char '\n'
  return $ BlastTabularHit (C.pack _queryId) (C.pack _subjectId) _seqIdentity _alignmentLength _misMatches _gapOpenScore _queryStart _queryEnd _hitSeqStart _hitSeqEnd _eValue _bitScore _subjectFrame (C.pack _querySeq) (C.pack _subjectSeq)
  
--IUPAC amino acid with gap
aminoacidLetters :: Char -> Bool
aminoacidLetters = inClass "ARNDCQEGHILMFPSTWYVBZX-"

--IUPAC nucleic acid characters with gap
nucleotideLetters :: Char -> Bool
nucleotideLetters = inClass "AGTCURYSWKMBDHVN-."

--IUPAC nucleic acid characters with gap
bioLetters :: Char -> Bool
bioLetters = inClass "ABCDEFGHIJKLMNOPQRSTVWXYZ.-"
