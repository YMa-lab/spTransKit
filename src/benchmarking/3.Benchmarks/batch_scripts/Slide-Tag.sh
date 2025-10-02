#!/bin/bash
#SBATCH -p gpu --gres=gpu:1
#SBATCH --partition=3090-gcondo
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=200G
#SBATCH -t 168:00:00
#SBATCH -o Slide-Tag.out
#SBATCH --job-name Slide-Tag

for id in Human_Cortex Human_Melanoma Human_Tonsil; do
    python3 Experiments/3.Benchmarks/benchmarking.py --sample_id $id --technology Slide-Tag --species "Homo sapiens" --num_hvg 500 --num_neighbors 10
done