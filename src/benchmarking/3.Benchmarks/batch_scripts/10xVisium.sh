#!/bin/bash
#SBATCH -p gpu
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=190G
#SBATCH -t 24:00:00
#SBATCH -o 10xVisium.out
#SBATCH --job-name 10xVisium

for id in DLPFC_151507 DLPFC_151508 DLPFC_151509 DLPFC_151510 DLPFC_151669 DLPFC_151670 DLPFC_151671 DLPFC_151672 DLPFC_151673 DLPFC_151674 DLPFC_151675 DLPFC_151676; do
    python3 Experiments/3.Benchmarks/benchmarking.py --sample_id $id --technology 10xVisium --species "Homo sapiens" --num_hvg 500 --num_neighbors 10
done