library(BayesSpace)
library(clue)
library(SingleCellExperiment)
library(ggplot2)
library(Seurat)
library(dplyr)
library(mclust)
library(Matrix)
library(zellkonverter)

args <- commandArgs(trailingOnly = TRUE)
sample.name <- as.numeric(args[1])
hvg_type <- as.character(args[2])
  
print(paste("Sample Name started:", sample.name))
dir.output <- file.path(paste0('/oscar/data/yma16/Project/spTransform/2.Downstream/results/10xVisium/DLPFC_hvg_', hvg_type, "_sdd_kmeans"), sample.name, 'bayesspace')

filepath_norm <- file.path(paste0('/oscar/data/yma16/Project/spTransform/0.NormalizedCounts/10xVisium/DLPFC_', sample.name), 'norm_counts')
norm_list <- list.files(filepath_norm, full.names = TRUE)
norm_list <- norm_list[file.info(norm_list)$isdir == FALSE]
norm_list <- basename(norm_list)
norm_list <- sub("\\.h5ad$", "", norm_list)
norm_list <- norm_list[! norm_list %in% c("scanpy_seurat", "scanpy_weinreb")]
norm_list <- c(norm_list, "original_bayesspace")

for (normalization in norm_list) {
  print(paste("Normalization Method started:", normalization))
  start_time <- Sys.time()
  dir.outpath <- file.path(dir.output, normalization)
  if(!dir.exists(file.path(dir.outpath))){
    dir.create(file.path(dir.outpath), recursive = TRUE)
  }
  if (!file.exists(file.path(dir.outpath, "10xVisium_DLPFC_run5_bayesspace_ari_score.txt"))){
    if (normalization == "original_bayesspace") {
      # load counts
      count_file <- file.path(paste0('/oscar/data/yma16/Project/spTransform/0.NormalizedCounts/10xVisium/DLPFC_', sample.name), 'norm_counts/raw.h5ad')
      counts <- readH5AD(count_file, use_hdf5 = TRUE)
      counts_matrix <- as(assay(counts, "X"), "CsparseMatrix")
      assay(counts, "X") <- counts_matrix
      order <- colnames(counts_matrix)
      rm(counts_matrix)
      gc()
      assay(counts, "logcounts") <- log1p(assay(counts, "X"))
      set.seed(101)
      dec <- scran::modelGeneVar(counts)
      top <- scran::getTopHVGs(dec, n = 2000)
      set.seed(102)
      counts <- scater::runPCA(counts, subset_row=top)
      rm(dec)
      rm(top)
      gc()      
      spatial_data <- as.Seurat(counts, counts = "X", data = "logcounts")
      rm(counts)
      gc()        
      }
    else{
      # load count
      count_file <- file.path(paste0('/oscar/data/yma16/Project/spTransform/0.NormalizedCounts/10xVisium/DLPFC_', sample.name), 'norm_counts',paste0(normalization, '.h5ad'))
      counts <- readH5AD(count_file, use_hdf5 = TRUE)
      counts_matrix <- as(assay(counts, "X"), "CsparseMatrix")

      # load hvg (cv or var) and then take the subset
      hvg_path <- file.path("/oscar/data/yma16/Project/spTransform/1.Evaluations/10xVisium", paste0("DLPFC_", sample.name), paste0(hvg_type, "_hvg_sdd"), paste0(normalization, ".csv"))

      hvg_genes_df <- read.csv(hvg_path, header = TRUE, stringsAsFactors = FALSE)
      hvg_genes <- hvg_genes_df[[2]]
      message(sprintf("Number of HVGs for sample %s, norm %s: %d",
                sample.name, normalization, length(hvg_genes)))
      counts_hvg <- counts_matrix[rownames(counts_matrix) %in% hvg_genes, , drop = FALSE]

      # load spatial
      spatial_data <- CreateSeuratObject(counts = counts_hvg, assay = "Spatial")
      order <- colnames(counts)
      rm(counts)
      rm(counts_matrix)
      rm(counts_hvg)
      gc() 
    }
    
    # load coords
    coord_file <- file.path(paste0('/oscar/data/yma16/Project/spTransform/0.NormalizedCounts/10xVisium/DLPFC_', sample.name), 'coordinates.csv')
    coord <- read.csv(coord_file, row.names = 1)
    spatial_data <- AddMetaData(spatial_data, metadata = coord$x, col.name = 'row')
    spatial_data <- AddMetaData(spatial_data, metadata = coord$y, col.name = 'col')
    rm(coord)
    gc()     
    
    spatial_dir <- file.path('/oscar/data/yma16/Project/spTransform/2.Downstream/datasets/10xVisium/DLPFC', sample.name, 'spatial')
    image <- Read10X_Image(image.dir = spatial_dir)
    
    # Attach the image
    spatial_data[["image"]] <- image
    if (normalization == "original_bayesspace") {
      DefaultAssay(spatial_data) <- "originalexp"
    }
    else {
      DefaultAssay(spatial_data) <- "Spatial"
    }
    rm(image)
    gc() 

    # Add annotations
    file_dic_gt <- '/oscar/data/yma16/Project/spTransform/2.Downstream/datasets/10xVisium/DLPFC/'
    df_anno <- read.csv(paste0(file_dic_gt, "layers.csv"), row.names = 1)
    df_anno_cp <- df_anno %>% filter(sample_id == as.numeric(sample.name))
    rm(df_anno)
    gc() 
    rownames(df_anno_cp) <- sub("\\..*$", "", rownames(df_anno_cp))
    df_anno_cp <- df_anno_cp[order, ]
    spatial_data <- AddMetaData(spatial_data, metadata = df_anno_cp$layer_guess_reordered, col.name = "ground_truth")
    ground_truth_values <- FetchData(spatial_data, vars = "ground_truth")
    ground_truth_values <- na.omit(ground_truth_values)
    spatial_data <- subset(spatial_data, cells = rownames(ground_truth_values))
    numCluster <- length(unique(spatial_data$ground_truth))
    rm(df_anno_cp)
    rm(ground_truth_values)
    gc() 
    
    # convert to a SingleCellExperiment object
    spatial_data_sce <- as.SingleCellExperiment(spatial_data)
    logcounts(spatial_data_sce) <- counts(spatial_data_sce)
    
    ## cluster
    for (i in 1:5) {
      tryCatch({
        set.seed(100+i)
        if (normalization == "original_bayesspace") {
          spatial_data_sce_cp <- spatial_data_sce
        }
        else {
          spatial_data_sce_cp <- scater::runPCA(spatial_data_sce,ncomponents = 15)
        }
        spatial_data_sce_cp <- spatialPreprocess(spatial_data_sce_cp, platform="Visium", skip.PCA=TRUE)
        sce <- spatialCluster(spatial_data_sce_cp, q=numCluster, d=15, platform='Visium',init.method="kmeans", model = "t", 
                                nrep=50000, gamma=3, save.chain=TRUE)
        comp_df <- as.data.frame(colData(sce))
        ari_score <- adjustedRandIndex(comp_df$ground_truth, comp_df$spatial.cluster)
        dir.outpath2 <- file.path(dir.outpath, paste0("10xVisium_DLPFC_run",i,"_bayesspace_"))
        write.table(comp_df,paste0(dir.outpath2, "domain.txt"))
        write(ari_score, file = paste0(dir.outpath2, "ari_score.txt"))
        rm(spatial_data_sce_cp)
        rm(sce)
        rm(comp_df)
        gc()
      },
      error = function(e) {
        message(sprintf("❌ Error for normalization '%s', run %d: %s",
                        normalization, i, e$message))
      })
    }
    rm(spatial_data)
    rm(image)
    rm(df_anno)
    rm(df_anno_cp)
    gc() 
    end_time <- Sys.time()
    execution_time <- end_time - start_time
    print(paste("Time taken:", execution_time))
  }
}
