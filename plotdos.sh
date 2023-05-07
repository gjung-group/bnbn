#!/bin/bash
for i in {-10..10}
do
cd ele_$i
sed -n 2p BLBL.mol> nn
num=`awk '{print $1}' nn`
tt=`echo "scale=0;$num/2" | bc -l`
bb=`echo "scale=0;$tt+1" | bc -l`
echo $num 
echo $tt
echo $bb

cat > plotdos.py <<EOF
import numpy as np
import matplotlib.pyplot as plt
import __tBG_dkl as tBG

LDOS = tBG.Load_object("./LDOS")
Egrid = tBG.Load_object("./E_LDOS")
top = 0
for i in range($tt):
    top += LDOS[i]
np.savetxt('energy',Egrid)
np.savetxt('top',top)

bot= 0
for j in range($bb,$num):
    bot += LDOS[j]
np.savetxt('bot',bot)
EOF
cd ../
done
