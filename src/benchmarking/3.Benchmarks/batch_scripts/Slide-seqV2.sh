#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p vnc
#SBATCH --mem=500G
#SBATCH -t 168:00:00
#SBATCH -o Slide-seqV2.out
#SBATCH --job-name Slide-seqV2

for id in Mouse_Cerebellum Mouse_Hippocampus; do
    python3 Experiments/3.Benchmarks/benchmarking.py --sample_id $id --technology Slide-seqV2 --species "Mus musculus" --num_hvg 500 --num_neighbors 10
done