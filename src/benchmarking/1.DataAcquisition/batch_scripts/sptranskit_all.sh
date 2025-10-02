#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=200G
#SBATCH -t 12:00:00
#SBATCH -o output.out
#SBATCH --job-name spTransKit

# 10X VISIUM

for id in DLPFC_151507 DLPFC_151508 DLPFC_151509 DLPFC_151510 DLPFC_151669 DLPFC_151670 DLPFC_151671 DLPFC_151672 DLPFC_151673 DLPFC_151674 DLPFC_151675 DLPFC_151676; do
    python3 Experiments/1.DataAcquisition/data_acquisition.py --sample_id $id --technology 10xVisium --counts data/10xVisium/${id}/counts.h5 --coordinates data/10xVisium/${id}/coordinates.csv
done

for id in Human_BCRA_CytAs Human_BCRA_FFPE Human_BCRA_VisiumFreshWT Human_BCRA_WTA Human_CESC_FFPE; do
    python3 Experiments/1.DataAcquisition/data_acquisition.py --sample_id $id --technology 10xVisium --counts data/10xVisium/${id}/counts.h5 --coordinates data/10xVisium/${id}/coordinates.csv 
done

for id in Human_COAD_11mmFFPE Human_COAD_CytAs Human_COAD_Intestine Human_COAD_WTA Human_GBM_WTA; do
    python3 Experiments/1.DataAcquisition/data_acquisition.py --sample_id $id --technology 10xVisium --counts data/10xVisium/${id}/counts.h5 --coordinates data/10xVisium/${id}/coordinates.csv
done

for id in Mouse_IPMN_15wksDox Mouse_IPMN_25wksDox; do
    python3 Experiments/1.DataAcquisition/data_acquisition.py --sample_id $id --technology 10xVisium --counts data/10xVisium/${id}/counts.h5 --coordinates data/10xVisium/${id}/coordinates.csv
done

for id in Human_IPMN_HG_1 Human_IPMN_HG_2 Human_IPMN_HG_3 Human_IPMN_LG_1 Human_IPMN_LG_2 Human_IPMN_LG_3 Human_IPMN_LG_4 Human_IPMN_LG_5 Human_IPMN_LG_6 Human_IPMN_LG_7 Human_IPMN_PDAC_1 Human_IPMN_PDAC_2 Human_IPMN_PDAC_3; do
    python3 Experiments/1.DataAcquisition/data_acquisition.py --sample_id $id --technology 10xVisium --counts data/10xVisium/${id}/counts.h5 --coordinates data/10xVisium/${id}/coordinates.csv
done

# SLIDE-SEQ

for id in Cerebellum_180430_1; do
    python3 src/benchmarking.py --sample_id $id --technology Slide-seq --counts data/Slide-seq/${id}/counts.csv --coordinates data/Slide-seq/${id}/coordinates.csv --species "Mus musculus" --genes data/Slide-seq/${id}/marker_genes.csv --pathways data/Slide-seq/${id}/marker_pathways.txt --num_hvg 500 --num_neighbors 10
done

for id in Human_Testis_Puck5 Human_Testis_Puck6; do
    python3 src/benchmarking.py --sample_id $id --technology Slide-seq --counts data/Slide-seq/${id}/counts.csv --coordinates data/Slide-seq/${id}/coordinates.csv --species "Homo sapiens" --genes data/Slide-seq/${id}/marker_genes.csv --pathways data/Slide-seq/${id}/marker_pathways.txt --num_hvg 500 --num_neighbors 10
done

for id in Mouse_Testis_Diabetes_1 Mouse_Testis_Diabetes_2 Mouse_Testis_Diabetes_3 Mouse_Testis_WT_1 Mouse_Testis_WT_2 Mouse_Testis_WT_3; do
    python3 src/benchmarking.py --sample_id $id --technology Slide-seq --counts data/Slide-seq/${id}/counts.csv --coordinates data/Slide-seq/${id}/coordinates.csv --species "Mus musculus" --genes data/Slide-seq/${id}/marker_genes.csv --pathways data/Slide-seq/${id}/marker_pathways.txt --num_hvg 500 --num_neighbors 10
done

# SLIDE-SEQV2

for id in Mouse_Cerebellum Mouse_Hippocampus; do
    python3 src/benchmarking.py --sample_id $id --technology Slide-seqV2 --counts data/Slide-seqV2/${id}/counts.csv --coordinates data/Slide-seqV2/${id}/coordinates.csv --species "Mus musculus" --genes data/Slide-seqV2/${id}/marker_genes.csv --pathways data/Slide-seqV2/${id}/marker_pathways.txt --num_hvg 500 --num_neighbors 10
done

# MERFISH

for id in MouseIleum; do
    python3 src/benchmarking.py --sample_id $id --technology Merfish --counts data/Merfish/${id}/counts.csv --coordinates data/Merfish/${id}/coordinates.csv --species "Mus musculus" --genes data/Merfish/${id}/marker_genes.csv --pathways data/Merfish/${id}/marker_pathways.txt --num_hvg 50 --num_neighbors 10
done

for id in HumanMTG_H18_250_1 HumanSTG_H19_250_1 HumanSTG_H20_250_1 HumanSTG_H20_250_2 HumanMTG_H22_250_1; do
    python3 src/benchmarking.py --sample_id $id --technology Merfish --counts data/Merfish/${id}/counts.csv --coordinates data/Merfish/${id}/coordinates.csv --species "Homo sapiens" --genes data/Merfish/${id}/marker_genes.csv --pathways data/Merfish/${id}/marker_pathways.txt --num_hvg 50 --num_neighbors 10
done

for id in HumanMTG_H18_4000_1 HumanMTG_H18_4000_2 HumanMTG_H18_4000_3 HumanSTG_H19_4000_1 HumanSTG_H19_4000_2 HumanSTG_H20_4000_1 HumanSTG_H20_4000_2 HumanSTG_H20_4000_3 HumanMTG_H22_4000_1 HumanMTG_H22_4000_2; do
    python3 src/benchmarking.py --sample_id $id --technology Merfish --counts data/Merfish/${id}/counts.csv --coordinates data/Merfish/${id}/coordinates.csv --species "Homo sapiens" --genes data/Merfish/${id}/marker_genes.csv --pathways data/Merfish/${id}/marker_pathways.txt --num_hvg 500 --num_neighbors 10
done

# SEQFISH

for id in MouseEmbryo; do
    python3 src/benchmarking.py --sample_id $id --technology seqFISH --counts data/seqFISH/${id}/counts.csv --coordinates data/seqFISH/${id}/coordinates.csv --species "Mus musculus" --genes data/seqFISH/${id}/marker_genes.csv --pathways data/seqFISH/${id}/marker_pathways.txt --num_hvg 500 --num_neighbors 10
done

# SEQFISH+

for id in Mouse_Olf_Bulb Mouse_SVZ; do
    python3 src/benchmarking.py --sample_id $id --technology seqFISH+ --counts data/seqFISH+/${id}/counts.csv --coordinates data/seqFISH+/${id}/coordinates.csv --species "Mus musculus" --genes data/seqFISH+/${id}/marker_genes.csv --pathways data/seqFISH+/${id}/marker_pathways.txt --num_hvg 500 --num_neighbors 10
done

# SLIDE-TAG

for id in Human_Cortex Human_Melanoma Human_Tonsil; do
    python3 src/benchmarking.py --sample_id $id --technology Slide-Tag --counts data/Slide-Tag/${id}/counts.csv --coordinates data/Slide-Tag/${id}/coordinates.csv --species "Homo sapiens" --genes data/Slide-Tag/${id}/marker_genes.csv --pathways data/Slide-Tag/${id}/marker_pathways.txt --num_hvg 500 --num_neighbors 10
done

# STEREO-SEQ

for id in Mouse_Brain_S; do
    python3 src/benchmarking.py --sample_id $id --technology Stereo-seq --counts data/Stereo-seq/${id}/counts.h5ad --coordinates data/Stereo-seq/${id}/coordinates.csv --species "Mus musculus" --genes data/Stereo-seq/${id}/marker_genes.csv --pathways data/Stereo-seq/${id}/marker_pathways.txt --num_hvg 500 --num_neighbors 10
done

# 10X XENIUM

for id in Human_LUAD Human_Brain Human_Breast_1 Human_Breast_2 Human_Skin; do
    python3 src/benchmarking.py --sample_id $id --technology 10x_Xenium --counts data/10x_Xenium/${id}/counts.h5 --coordinates data/10x_Xenium/${id}/coordinates.csv --species "Homo sapiens" --genes data/10x_Xenium/${id}/marker_genes.csv --pathways data/10x_Xenium/${id}/marker_pathways.txt --num_hvg 50 --num_neighbors 10
done

for id in Mouse_Brain; do
    python3 src/benchmarking.py --sample_id $id --technology 10x_Xenium --counts data/10x_Xenium/${id}/counts.h5 --coordinates data/10x_Xenium/${id}/coordinates.csv --species "Mus musculus" --genes data/10x_Xenium/${id}/marker_genes.csv --pathways data/10x_Xenium/${id}/marker_pathways.txt --num_hvg 50 --num_neighbors 10
done

# 10 XENIUM PRIME

for id in Human_Prostate Human_Skin_P; do
    python3 src/benchmarking.py --sample_id $id --technology 10x_Xenium_Prime --counts data/10x_Xenium_Prime/${id}/counts.h5 --coordinates data/10x_Xenium_Prime/${id}/coordinates.csv --species "Homo sapiens" --genes data/10x_Xenium_Prime/${id}/marker_genes.csv --pathways data/10x_Xenium_Prime/${id}/marker_pathways.txt --num_hvg 500 --num_neighbors 10
done

# NANOSTRING COSMX

for id in Human_NSCLC_5_1 Human_NSCLC_5_2 Human_NSCLC_5_3 Human_NSCLC_6 Human_NSCLC_9_1 Human_NSCLC_9_2 Human_NSCLC_12 Human_NSCLC_13; do
    python3 src/benchmarking.py --sample_id $id --technology NanoString_CosMx --counts data/NanoString_CosMx/${id}/counts.csv --coordinates data/NanoString_CosMx/${id}/coordinates.csv --species "Homo sapiens" --genes data/NanoString_CosMx/${id}/marker_genes.csv --pathways data/NanoString_CosMx/${id}/marker_pathways.txt --num_hvg 100 --num_neighbors 10
done

# OPEN-ST

for id in Mouse_Embryo_1 Mouse_Embryo_2; do
    python3 src/benchmarking.py --sample_id $id --technology Open-ST --counts data/Open-ST/${id}/counts.h5ad --coordinates data/Open-ST/${id}/coordinates.csv --species "Mus musculus" --genes data/Open-ST/${id}/marker_genes.csv --pathways data/Open-ST/${id}/marker_pathways.txt --num_hvg 500 --num_neighbors 10
done

for id in Human_HNSCC Human_MLN_2 Human_MLN_3 Human_MLN_4 Human_MLN_5; do
    python3 src/benchmarking.py --sample_id $id --technology Open-ST --counts data/Open-ST/${id}/counts.h5ad --coordinates data/Open-ST/${id}/coordinates.csv --species "Homo sapiens" --genes data/Open-ST/${id}/marker_genes.csv --pathways data/Open-ST/${id}/marker_pathways.txt --num_hvg 500 --num_neighbors 10
done

for id in Human_MLN_6 Human_MLN_7 Human_MLN_9 Human_MLN_11 Human_MLN_18; do
    python3 src/benchmarking.py --sample_id $id --technology Open-ST --counts data/Open-ST/${id}/counts.h5ad --coordinates data/Open-ST/${id}/coordinates.csv --species "Homo sapiens" --genes data/Open-ST/${id}/marker_genes.csv --pathways data/Open-ST/${id}/marker_pathways.txt --num_hvg 500 --num_neighbors 10
done

for id in Human_MLN_19 Human_MLN_23 Human_MLN_24 Human_MLN_25; do
    python3 src/benchmarking.py --sample_id $id --technology Open-ST --counts data/Open-ST/${id}/counts.h5ad --coordinates data/Open-ST/${id}/coordinates.csv --species "Homo sapiens" --genes data/Open-ST/${id}/marker_genes.csv --pathways data/Open-ST/${id}/marker_pathways.txt --num_hvg 500 --num_neighbors 10
done

for id in Human_MLN_28 Human_MLN_33 Human_MLN_34 Human_MLN_36; do
    python3 src/benchmarking.py --sample_id $id --technology Open-ST --counts data/Open-ST/${id}/counts.h5ad --coordinates data/Open-ST/${id}/coordinates.csv --species "Homo sapiens" --genes data/Open-ST/${id}/marker_genes.csv --pathways data/Open-ST/${id}/marker_pathways.txt --num_hvg 500 --num_neighbors 10
done