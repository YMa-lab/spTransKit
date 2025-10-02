#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p viz
#SBATCH --mem=700G
#SBATCH -t 168:00:00
#SBATCH -o 10x_Xenium_Prime.out
#SBATCH --job-name 10x_Xenium_Prime

for id in Human_Skin_P; do
    python3 Experiments/3.Benchmarks/benchmarking.py --sample_id $id --technology 10x_Xenium_Prime --species "Homo sapiens" --num_hvg 500 --num_neighbors 10
done