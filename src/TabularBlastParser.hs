{-# LANGUAGE OverloadedStrings #-}

-- This attoparsec module
module TabularBlastParser (module TabularBlastParser)
where
import Data.Attoparsec.ByteString.Char8
import Data.Word
import qualified Data.ByteString.Char8 as C
import qualified Data.Vector as V

data BlastTabularResult = BlastTabularResult
  { blastProgram :: C.ByteString,
    blastQuery :: C.ByteString,
    blastDatabase :: C.ByteString,
    blastHitNumber :: Int,
    hitLines :: V.Vector BlastTabularHit
  }
  deriving (Show, Eq)

parseTabularBlastHeader :: Parser BlastTabularResult
parseTabularBlastHeader = do
  string "# "
  _blastProgram <-  many1 (notChar '\n')
  string "\n# Query: "
  _blastQuery <- many1 (notChar '\n')
  string "\n# Database: "
  _blastDatabase <- many1 (notChar '\n')
  string "\n# "
  _blastHitNumber <- decimal
  string " hits found\n"
  _tabularHit <- many1 parseBlastTabularHit
  return $ BlastTabularResult (C.pack _blastProgram) (C.pack _blastQuery) (C.pack _blastDatabase) _blastHitNumber (V.fromList _tabularHit)

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

parseBlastTabularHit :: Parser BlastTabularHit
parseBlastTabularHit = do
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
