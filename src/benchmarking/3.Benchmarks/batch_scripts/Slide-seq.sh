#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p vnc
#SBATCH --mem=700G
#SBATCH -t 168:00:00
#SBATCH -o Slide-seq.out
#SBATCH --job-name Slide-seq

for id in Cerebellum_180430_1 Human_Testis_Puck5 Human_Testis_Puck6; do
    python3 Experiments/3.Benchmarks/benchmarking.py --sample_id $id --technology Slide-seq --species "Homo sapiens" --num_hvg 500 --num_neighbors 10
done

for id in Mouse_Testis_Diabetes_1 Mouse_Testis_Diabetes_2 Mouse_Testis_Diabetes_3 Mouse_Testis_WT_1 Mouse_Testis_WT_2 Mouse_Testis_WT_3; do
    python3 Experiments/3.Benchmarks/benchmarking.py --sample_id $id --technology Slide-seq --species "Mus musculus" --num_hvg 500 --num_neighbors 10
done