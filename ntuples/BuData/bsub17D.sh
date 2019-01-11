#!/bin/bash
i=0;
max=35;
n=1000000;
while [ "$i" -le "$max" ]; do
  mkdir "s17D_$i";
  cd "s17D_$i";
  skip=$(python -c "print 0+$n*$i");
  echo $'#!/bin/sh' > script.sh
  echo $'#BSUB -o test.log' >> script.sh
  echo $'eval `scram runtime -sh`' >> script.sh
  echo "pdTreeAnalyze /lustre/cmswork/abragagn/ntuList/charmonium2017Lists/Charmonium_Run2017D-31Mar2018-v1_MINIAOD_DCAP.list hist$i.root -v outputFile ntu$i.root -v histoMode RECREATE -v use_gen f -v mvaMethod BDTMuonID2017woIPwIso -v muonIdWpBarrel 0.8910 -v muonIdWpEndcap 0.8925 -v useHLT t -v process BuJPsiK -n $n -s $skip" >> script.sh
  echo "" >> script.sh
  bsub < script.sh;
  cd ..;
  i=`expr "$i" + 1`;
done