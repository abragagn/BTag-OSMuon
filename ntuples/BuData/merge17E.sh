#!/bin/bash
i=0;
max=60;
echo $'#!/bin/sh' > park.sh
echo -n "hadd ntuBuData2017E.root" >> park.sh
while [ "$i" -le "$max" ]; do
  echo -n " s17E_$i/ntu$i" >> park.sh
  echo -n ".root" >> park.sh
  i=`expr "$i" + 1`;
done
echo " " >> park.sh
bash park.sh;
rm park.sh;
mv ntuBuData2017E.root ../
