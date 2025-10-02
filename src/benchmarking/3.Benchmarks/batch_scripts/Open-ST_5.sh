#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p bigmem
#SBATCH --mem=700G
#SBATCH -t 48:00:00
#SBATCH -o Open-ST_5.out
#SBATCH --job-name Open-ST_5

for id in Human_MLN_36; do
    python3 Experiments/3.Benchmarks/benchmarking.py --sample_id $id --technology Open-ST --species "Homo sapiens" --num_hvg 500 --num_neighbors 10
done