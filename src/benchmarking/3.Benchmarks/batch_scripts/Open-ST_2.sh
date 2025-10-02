#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p bigmem
#SBATCH --mem=700G
#SBATCH -t 48:00:00
#SBATCH -o Open-ST_2_2.out
#SBATCH --job-name Open-ST_2_2

for id in Human_MLN_2; do
    python3 Experiments/3.Benchmarks/benchmarking.py --sample_id $id --technology Open-ST --species "Homo sapiens" --num_hvg 500 --num_neighbors 10
done