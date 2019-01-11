#!/bin/sh
#BSUB -o test231.log
eval `scram runtime -sh`
python trainingOsMuon.py 231 100 3 0.25
