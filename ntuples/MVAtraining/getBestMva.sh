#!/bin/bash
i=1
j=1
k=0
max_i=5
max_j=4
max_k=2
while [ "$i" -le "$max_i" ]; do
    while [ "$j" -le "$max_j" ]; do
        while [ "$k" -le "$max_k" ]; do
            mkdir "mva$i$j$k";
            cd "mva$i$j$k";
            echo $'#!/bin/sh' > script.sh
            echo $'#BSUB -o test$i$j$ki.log' >> script.sh
            echo "root -l -b -q '../fitMVA.C+(\"ntuBsMC2017.root\", \"DNNOsMuonHLTJpsiMu_test$i$j$ki\", ,\"CREATE\", false, false, -1)'" >> script.sh
            echo "" >> script.sh
            bsub < script.sh;
            cd ..;
            k=`expr "$k" + 1`
        done
        k=0
        j=`expr "$j" + 1`
    done
    j=1
    i=`expr "$i" + 1`
done
