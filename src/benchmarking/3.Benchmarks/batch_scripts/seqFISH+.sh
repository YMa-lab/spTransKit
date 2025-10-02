#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p batch
#SBATCH --mem=200G
#SBATCH -t 168:00:00
#SBATCH -o seqFISH+.out
#SBATCH --job-name seqFISH+

# seqFISH+

for id in Mouse_Olf_Bulb Mouse_SVZ; do
    python3 Experiments/3.Benchmarks/benchmarking.py --sample_id $id --technology seqFISH+ --species "Mus musculus" --num_hvg 500 --num_neighbors 10
done