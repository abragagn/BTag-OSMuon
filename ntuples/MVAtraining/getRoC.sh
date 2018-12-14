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
            stuff=$(tail -n 17 testLayers$i$j$k.log | head -n 1 | tail -c 11)
            echo $stuff >> ROC.txt
            k=`expr "$k" + 1`
        done
        k=0
        j=`expr "$j" + 1`
    done
    j=1
    i=`expr "$i" + 1`
done
