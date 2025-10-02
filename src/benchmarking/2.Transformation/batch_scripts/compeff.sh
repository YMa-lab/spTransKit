#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p vnc
#SBATCH --mem=700G
#SBATCH -t 96:00:00
#SBATCH -o compeff.out
#SBATCH --job-name compeff

module load r/4.4
module load imagemagick/7.1.1-3-ex4k4u2

Rscript Experiments/2.Transformation/compeff_r.R