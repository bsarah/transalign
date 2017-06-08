module Main where

import System.Environment (getArgs)

import Transalign.BlastCache

main :: IO ()
main = mapM_ buildCache =<< getArgs
