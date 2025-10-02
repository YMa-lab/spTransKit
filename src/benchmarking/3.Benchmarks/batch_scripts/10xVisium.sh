#!/bin/bash
#SBATCH -p gpu
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=190G
#SBATCH -t 24:00:00
#SBATCH -o 10xVisium.out
#SBATCH --job-name 10xVisium

for id in Human_IPMN_LG_6; do
    python3 Experiments/3.Benchmarks/benchmarking.py --sample_id $id --technology 10xVisium --species "Homo sapiens" --num_hvg 500 --num_neighbors 10
done