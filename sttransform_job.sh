#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=1000G
#SBATCH -t 24:00:00
#SBATCH -o output.out
#SBATCH --job-name stTransform

python3 src/transformations.py --counts data/151673/matrix.csv --coordinates data/151673/coordinates.txt --genes data/151673/marker_genes.csv --pathways data/151673/marker_pathways.txt --sample_id 151673
python3 src/transformations.py --counts data/151674/matrix.csv --coordinates data/151674/coordinates.txt --genes data/151674/marker_genes.csv --pathways data/151674/marker_pathways.txt --sample_id 151674
python3 src/transformations.py --counts data/151675/matrix.csv --coordinates data/151675/coordinates.txt --genes data/151675/marker_genes.csv --pathways data/151675/marker_pathways.txt --sample_id 151675
python3 src/transformations.py --counts data/151676/matrix.csv --coordinates data/151676/coordinates.txt --genes data/151676/marker_genes.csv --pathways data/151676/marker_pathways.txt --sample_id 151676