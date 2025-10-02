#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p gpu
#SBATCH --mem=190G
#SBATCH -t 96:00:00
#SBATCH -o 10xVisium.out
#SBATCH --job-name 10xVisium

module load r/4.4
module load imagemagick/7.1.1-3-ex4k4u2

for id in DLPFC_151507 DLPFC_151508 DLPFC_151509 DLPFC_151510 DLPFC_151669 DLPFC_151670 DLPFC_151671 DLPFC_151672 DLPFC_151673 DLPFC_151674 DLPFC_151675 DLPFC_151676; do
    python3 Experiments/2.Transformation/transformations.py --sample_id $id --technology 10xVisium --num_hvg 500
    Rscript Experiments/2.Transformation/transformations_in_r.R --sample_id $id --technology 10xVisium
    Rscript Experiments/2.Transformation/r_trans_subcounts.R --sample_id $id --technology 10xVisium
    python3 Experiments/2.Transformation/integrate_data.py --sample_id $id --technology 10xVisium --num_hvg 500
done