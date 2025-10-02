#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p bigmem
#SBATCH --mem=700G
#SBATCH -t 48:00:00
#SBATCH -o NanoString_CosMx.out
#SBATCH --job-name NanoString_CosMx

for id in Human_NSCLC_9_2; do
    python3 Experiments/3.Benchmarks/benchmarking.py --sample_id $id --technology NanoString_CosMx --species "Homo sapiens" --num_hvg 100 --num_neighbors 10
done