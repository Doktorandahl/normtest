#!/bin/bash -l
#SBATCH -A naiss2023-5-551
#SBATCH -J normtestbw
#SBATCH -p main
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -o normaltest1bw.out
Rscript normal_test_mv_dardel_bw.R
