#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p batch
#SBATCH --mem=490G
#SBATCH -t 168:00:00
#SBATCH -o Open-ST_1.out
#SBATCH --job-name Open-ST_1

for id in Mouse_Embryo_1 Mouse_Embryo_2; do
    python3 Experiments/3.Benchmarks/benchmarking.py --sample_id $id --technology Open-ST --species "Mus musculus" --num_hvg 500 --num_neighbors 10
done

for id in Human_HNSCC; do
    python3 Experiments/3.Benchmarks/benchmarking.py --sample_id $id --technology Open-ST --species "Homo sapiens" --num_hvg 500 --num_neighbors 10
done