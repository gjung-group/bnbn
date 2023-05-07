#!/bin/bash
for i in {1..30}
do
cp atom-pos.py topbot.sh mm_$i
cd mm_$i
./topbot.sh
python atom-pos.py
cd ../
done
