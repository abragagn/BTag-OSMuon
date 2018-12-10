#!/bin/bash
i=0;
max=20;
n=1000000;
while [ "$i" -le "$max" ]; do
  mkdir "s16_$i";
  cd "s16_$i";
  skip=$(python -c "print 0+$n*$i");
  echo $'#!/bin/sh' > script.sh
  echo $'#BSUB -o test.log' >> script.sh
  echo $'eval `scram runtime -sh`' >> script.sh
  echo "pdTreeAnalyze /lustre/cmswork/abragagn/ntuList/MC2016Lists/BuToJpsiK_BMuonFilter_2016_DCAP.list hist$i.root -v outputFile ntu$i.root -v histoMode RECREATE -v use_gen t -v mvaMethod DNNGlobal2016woIPwIso -v muonIdWpBarrel 0.26 -v muonIdWpEndcap 0.60 -v useHLT t -v process BuJPsiK -n $n -s $skip" >> script.sh
  echo "" >> script.sh
  bsub < script.sh;
  cd ..;
  i=`expr "$i" + 1`;
done
