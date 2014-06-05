import BlastCache
import System.Environment (getArgs)

main :: IO ()
main = mapM_ buildCache =<< getArgs
