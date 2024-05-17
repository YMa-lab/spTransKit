#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=128G
#SBATCH -t 2:00:00
#SBATCH -o output.out
#SBATCH --job-name stTransform

python3 src/transformations.py --counts data/151673/matrix.csv --coordinates data/151673/coordinates.txt --genes data/151673/marker_genes.csv --pathways data/151673/marker_pathways.txt