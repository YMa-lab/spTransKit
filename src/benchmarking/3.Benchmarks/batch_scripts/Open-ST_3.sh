#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p vnc
#SBATCH --mem=900G
#SBATCH -t 168:00:00
#SBATCH -o Open-ST_3_3.out
#SBATCH --job-name Open-ST_3_3

for id in Human_MLN_25; do
    python3 Experiments/3.Benchmarks/benchmarking.py --sample_id $id --technology Open-ST --species "Homo sapiens" --num_hvg 500 --num_neighbors 10
done