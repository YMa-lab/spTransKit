#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p viz
#SBATCH --mem=800G
#SBATCH -t 168:00:00
#SBATCH -o Open-ST.out
#SBATCH --job-name Open-ST

for id in Human_MLN_11 Human_MLN_25; do
    python3 Experiments/3.Benchmarks/benchmarking.py --sample_id $id --technology Open-ST --species "Homo sapiens" --num_hvg 500 --num_neighbors 10
done