#!/bin/bash -l
#SBATCH -A naiss2023-5-551
#SBATCH -J normtest1
#SBATCH -p main
#SBATCH -t 00:10:00
#SBATCH -N 1
#SBATCH -o normaltest1.out
Rscript normal_test_mv_dardel.R
