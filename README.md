<img src="logo/sptranskit.png" width="75px" align="left"/> spTransKit: A Toolkit for Transformation Methods in Spatial Transcriptomics 
============

This repository provides the code for the 16 transformation methods evaluated in the study: A Comprehensive Benchmarking and Practical Guide to Transformation Methods for Spatial Transcriptomics and Downstream Analyses. The methods are designed to be easily called within any spatial transcriptomics analysis pipeline.

# Table of Contents
- [Background](#Background)
- [Transformations](#Transformations)

# Background

Spatial resolved transcriptomics (SRT) allows for the localization of gene expression to specific regions of tissue, aiding in the investigation of spatially dependent biological phenomena. Due to the many advantages of SRT over other transcriptomics technologies, several computational methods have been designed to analyze spatial transcriptomics data and extract biologically relevant spatial information. Despite the diversity of these methods, all pipelines typically begin with preprocessing of the raw expression data. Preprocessing is required to correct for the technical noise introduced by the spatial transcriptomics platform, which often obscures underlying biological signals.

# Transformations
| Name | Category | Function | Description |
|   :---:   |   :---:   |   :---:   |   :---:   |
| y/s | Library Size Factor-Based | size | Adjusts gene counts by the library size factor for each spatial location. |
| CPM | Library Size Factor-Based | cpm | Adjusts gene counts by the counts per million (CPM) library size factor for each spatial location. |
| scanpy Weinreb | Library Size Factor-Based | weinreb | Adjusts gene counts by the library size factor for each spatial location. |
| scanpy Zheng | Library Size Factor-Based | zheng | Adjusts gene counts by the library size factor for each spatial location. |
| TMM | Library Size Factor-Based | tmm | Adjusts gene counts by the library size factor for each spatial location. |
| DESeq2 | Library Size Factor-Based | deseq2 | Adjusts gene counts by the library size factor for each spatial location. |
