#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p gpu
#SBATCH --mem=190G
#SBATCH -t 96:00:00
#SBATCH -o Stereo-seq.out
#SBATCH --job-name Stereo-seq

module load r/4.4
module load imagemagick/7.1.1-3-ex4k4u2

for id in Mouse_Kidney Mouse_Lung Mouse_Testis Mouse_Tongue; do
    python3 Experiments/2.Transformation/transformations.py --sample_id $id --technology Stereo-seq --num_hvg 500
    Rscript Experiments/2.Transformation/transformations_in_r.R --sample_id $id --technology Stereo-seq
    Rscript Experiments/2.Transformation/r_trans_subcounts.R --sample_id $id --technology Stereo-seq
    python3 Experiments/2.Transformation/integrate_data.py --sample_id $id --technology Stereo-seq --num_hvg 500
done