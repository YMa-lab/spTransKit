#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p batch
#SBATCH --mem=300G
#SBATCH -t 96:00:00
#SBATCH -o subcounts.out
#SBATCH --job-name subcounts

python3 Experiments/3.Benchmarks/test.py