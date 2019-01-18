#!/bin/bash
i=1
j=1
k=0
max_i=5
max_j=4
max_k=2
rm bestMVA.txt;
while [ "$i" -le "$max_i" ]; do
    while [ "$j" -le "$max_j" ]; do
        while [ "$k" -le "$max_k" ]; do
            cd "mva$i$j$k";
            stuff=$(tail -n 8 test.log | head -n 1 | tail -c 9)
            echo -n "$i$j$k " >> ../bestMVA.txt
            echo $stuff >> ../bestMVA.txt
            cd ..;
            k=`expr "$k" + 1`
        done
        k=0
        j=`expr "$j" + 1`
    done
    j=1
    i=`expr "$i" + 1`
done

cat bestMVA.txt;
