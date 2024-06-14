#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=200G
#SBATCH -t 12:00:00
#SBATCH -o output.out
#SBATCH --job-name stTransform

for id in 151507 151508 151509 151510 151669 151670 151671 151672 151673 151674 151675 151676; do
    python3 src/transformations.py --sample_id $id --technology 10xVisium --counts data/10xVisium/${id}/counts.h5 --coordinates data/10xVisium/${id}/coordinates.csv  --genes data/10xVisium/${id}/marker_genes.csv --pathways data/10xVisium/${id}/marker_pathways.txt --num_hvg 500
done

for id in 180430_1; do
    python3 src/transformations.py --sample_id $id --technology Slide-seq --counts data/Slide-seq/${id}/counts.csv --coordinates data/Slide-seq/${id}/coordinates.csv  --genes data/Slide-seq/${id}/marker_genes.csv --pathways data/Slide-seq/${id}/marker_pathways.txt --num_hvg 500
done

for id in MouseHypothalamus MouseIleum HumanMTG_H18 HumanSTG_H19 HumanSTG_H20 HumanMTG_H22; do
    python3 src/transformations.py --sample_id $id --technology Merfish --counts data/Merfish/${id}/counts.csv --coordinates data/Merfish/${id}/coordinates.csv  --genes data/Merfish/${id}/marker_genes.csv --pathways data/Merfish/${id}/marker_pathways.txt --num_hvg 50
done

for id in MouseOlfactoryBulb; do
    python3 src/transformations.py --sample_id $id --technology Stereo-seq --counts data/Stereo-seq/${id}/counts.h5ad --coordinates data/Stereo-seq/${id}/coordinates.h5  --genes data/Stereo-seq/${id}/marker_genes.csv --pathways data/Stereo-seq/${id}/marker_pathways.txt --num_hvg 500
done