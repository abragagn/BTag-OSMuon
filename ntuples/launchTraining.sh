#!/bin/sh
#BSUB -o test.log
eval `scram runtime -sh`
python trainingOsMuon.py 2016 5
