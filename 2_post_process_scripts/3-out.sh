#!/bin/bash
for i in {1..30}
do
cp getOutput.py mm_$i
cd mm_$i
echo "out" $i
python getOutput.py 
cd ../
done
