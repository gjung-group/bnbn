#!/bin/bash 
#SBATCH -J doaa31405ele
#SBATCH -o results.o%j
#SBATCH -p mem740
#SBATCH -N 1
#SBATCH -n 32
##SBATCH -w n012

date
mpirun -np 32  python3 mkdos.py
date
