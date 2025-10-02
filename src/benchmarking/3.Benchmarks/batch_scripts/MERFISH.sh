#!/bin/bash
#SBATCH -p batch
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=200G
#SBATCH -t 168:00:00
#SBATCH -o MERFISH.out
#SBATCH --job-name MERFISH

for id in HumanMTG_H22_4000_2; do
    python3 Experiments/3.Benchmarks/benchmarking.py --sample_id $id --technology Merfish --species "Homo sapiens" --num_hvg 500 --num_neighbors 10
done