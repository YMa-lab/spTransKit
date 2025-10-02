#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p gpu
#SBATCH --mem=190G
#SBATCH -t 48:00:00
#SBATCH -o Stereo-seq.out
#SBATCH --job-name Stereo-seq

for id in Mouse_Tongue; do
    python3 Experiments/3.Benchmarks/benchmarking.py --sample_id $id --technology Stereo-seq --species "Mus musculus" --num_hvg 500 --num_neighbors 10
done