#!/bin/bash
i=1
j=1
k=0
max_i=5
max_j=4
max_k=2
while [ "$i" -le "$max_i" ]; do
    while [ "$j" -le "$max_j" ]; do
        size=$(python -c "print 50*$j")
        while [ "$k" -le "$max_k" ]; do
            drop=$(python -c "print 0.25*$k")
            echo "#!/bin/sh" > script.sh
            echo "#BSUB -o testLayers$i$j$k.log" >> script.sh
            echo "python trainingOsMuon.py $i$j$k $size $i $drop" >> script.sh
            bsub < script.sh;
            rm script.sh;
            k=`expr "$k" + 1`
        done
        k=0
        j=`expr "$j" + 1`
    done
    j=1
    i=`expr "$i" + 1`
done
