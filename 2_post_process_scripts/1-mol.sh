#!/bin/bash
for i in {1..30}
do
cp 1-fracToCart.py 2-charg.py atom-pos.py BN.drip BN.extep lammps.in lam.sh mm_$i
cp hei.sh 3-rep.sh mm_$i
cd mm_$i
python 1-fracToCart.py BNBN.frac
./3-rep.sh
python 2-charg.py 
cd ../
done
