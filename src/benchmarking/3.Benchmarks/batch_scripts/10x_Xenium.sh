#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p vnc
#SBATCH --mem=700G
#SBATCH -t 168:00:00
#SBATCH -o 10x_Xenium.out
#SBATCH --job-name 10x_Xenium

for id in Human_Skin; do
    python3 Experiments/3.Benchmarks/benchmarking.py --sample_id $id --technology 10x_Xenium --species "Homo sapiens" --num_hvg 50 --num_neighbors 10
done

for id in Mouse_Brain; do
    python3 Experiments/3.Benchmarks/benchmarking.py --sample_id $id --technology 10x_Xenium --species "Mus musculus" --num_hvg 50 --num_neighbors 10
done