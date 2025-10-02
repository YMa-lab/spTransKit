library(reticulate)
use_condaenv("giotto_env", conda = "/users/yzhu194/.local/share/r-miniconda/bin/conda", required = TRUE)
reticulate::py_config()
library(Giotto)
library(Matrix)
library(SingleCellExperiment)
anndata_py <- import("anndata", convert = FALSE)

technology <- 'seqFISH+'
sample.name <- 'Mouse_SVZ'
sample.name2 <- 'Mouse_SVZ'
hvg_type <- 'july_set'

print(paste("Technology started:", technology))
print(paste("Sample Name started:", sample.name))

dir.output <- file.path('/oscar/data/yma16/Project/spTransform/2.Downstream/results/seqFISH+/', technology, sample.name, paste0('Giotto_SCG_',hvg_type))
coord_file <- file.path('/oscar/data/yma16/Project/spTransform/0.NormalizedCounts', technology, sample.name2, 'coordinates.csv')

coord <- read.csv(coord_file, row.names = 1)
filepath_norm <- file.path('/oscar/data/yma16/Project/spTransform/0.NormalizedCounts/', technology, sample.name, 'norm_counts')
norm_list <- list.files(filepath_norm, full.names = TRUE)
norm_list <- norm_list[file.info(norm_list)$isdir == FALSE]
norm_list <- basename(norm_list)
norm_list <- sub("\\.h5ad$", "", norm_list)
norm_list <- c(norm_list, "original_giotto")
norm_list <- norm_list[!grepl("raw\\.csv$", basename(norm_list))]

final_spatial_genes_df <- read.csv("/oscar/data/yma16/Project/spTransform/2.Downstream/results/summary/SCG/all_spatial_genes_new.csv")
final_spatial_genes <- unique(unlist(final_spatial_genes_df))
norm_list <- c('original_giotto')
for (normalization in norm_list) {
  print(paste("Processing with final gene set for normalization:", normalization))
  dir.outpath <- file.path(dir.output, normalization)
  
  # Load expression matrix
  if (normalization == "original_giotto") {
    count_file <- file.path(filepath_norm, 'raw.h5ad')
    adata <- anndata_py$read_h5ad(count_file, backed = FALSE)
    X_mat <- py_to_r(adata$X)
    expr_matrix <- Matrix::Matrix(t(X_mat), sparse = TRUE)
    gene_names <- unlist(py_to_r(adata$var_names$tolist()))
    cell_names <- unlist(py_to_r(adata$obs_names$tolist()))
    rownames(expr_matrix) <- gene_names
    colnames(expr_matrix) <- cell_names
    expr_dense <- as.matrix(expr_matrix)
  }
  else {
    # load count
    count_file <- file.path(filepath_norm, paste0(normalization, '.h5ad'))
    adata <- anndata_py$read_h5ad(count_file, backed = FALSE)
    X_mat <- py_to_r(adata$X)
    expr_matrix <- Matrix::Matrix(t(X_mat), sparse = TRUE)
    gene_names <- unlist(py_to_r(adata$var_names$tolist()))
    cell_names <- unlist(py_to_r(adata$obs_names$tolist()))
    rownames(expr_matrix) <- gene_names
    colnames(expr_matrix) <- cell_names
    expr_dense <- as.matrix(expr_matrix)
  }
  
  instrs = createGiottoInstructions(save_plot = TRUE, 
                                    show_plot = FALSE,
                                    save_dir = dir.outpath)
  my_giotto_object = createGiottoObject(raw_exprs = expr_dense,
                                        spatial_locs = coord,
                                        instructions = instrs)
  
  # Filter and normalize as needed
  if (normalization == "original_giotto") {
    my_giotto_object <- filterGiotto(gobject = my_giotto_object, 
                                     expression_threshold = 0.5, 
                                     gene_det_in_min_cells = 20, 
                                     min_det_genes_per_cell = 0)
    my_giotto_object <- normalizeGiotto(gobject = my_giotto_object)
  } else {
    my_giotto_object <- filterGiotto(gobject = my_giotto_object, 
                                     expression_threshold = 0, 
                                     gene_det_in_min_cells = 0, 
                                     min_det_genes_per_cell = 0)
    my_giotto_object@norm_expr = expr_dense
  }
  
  # Create spatial network
  my_giotto_object = createSpatialNetwork(gobject = my_giotto_object, minimum_k = 2)
  
  ## run spatial gene co-expression method
  # 1. calculate spatial correlation scores 
  spat_cor_netw_DT = detectSpatialCorGenes(my_giotto_object,
                                           method = 'network', 
                                           spatial_network_name = 'Delaunay_network',
                                           subset_genes = final_spatial_genes,
                                           cor_method = "pearson")
  
  # 2. cluster correlation scores
  spat_cor_netw_DT = clusterSpatialCorGenes(spat_cor_netw_DT, 
                                            name = 'spat_netw_clus', k = 8)
  heatmSpatialCorGenes(my_giotto_object, spatCorObject = spat_cor_netw_DT, 
                       use_clus_name = 'spat_netw_clus')
  
  # rank spatial correlation clusters based on how similar they are
  netw_ranks = rankSpatialCorGroups(my_giotto_object, 
                                    spatCorObject = spat_cor_netw_DT, 
                                    use_clus_name = 'spat_netw_clus')
  
  # extract information about clusters
  top_netw_spat_cluster = showSpatialCorGenes(spat_cor_netw_DT, 
                                              use_clus_name = 'spat_netw_clus',
                                              selected_clusters = 6, 
                                              show_top_genes = 1)
  
  cluster_genes_DT = showSpatialCorGenes(spat_cor_netw_DT, 
                                         use_clus_name = 'spat_netw_clus',
                                         show_top_genes = 1)
  
  cluster_genes = cluster_genes_DT$clus; names(cluster_genes) = cluster_genes_DT$gene_ID
  write.csv(cluster_genes, file.path(dir.outpath, "giotto_cluster_genes.csv"))
  
  # create spatial metagenes and visualize
  my_giotto_object = createMetagenes(my_giotto_object, gene_clusters = cluster_genes, name = 'cluster_metagene')
  spatCellPlot(my_giotto_object,
               spat_enr_names = 'cluster_metagene',
               cell_annotation_values = netw_ranks$clusters,
               point_size = 1.5, cow_n_col = 3)
  
  gene_order <- spat_cor_netw_DT$gene_order
  mb_spark_est = spat_cor_netw_DT$cor_DT
  mb_spark_est$row = factor(mb_spark_est$gene_ID,levels = gene_order) %>%
    as.numeric()
  mb_spark_est$col = factor(mb_spark_est$variable,levels = gene_order) %>%
    as.numeric()
  
  
  mb_spark_corr = matrix(0,nrow = length(gene_order),ncol = length(gene_order))
  mb_spark_corr[as.matrix(mb_spark_est[,c("row","col")])] = mb_spark_est$spat_cor
  rownames(mb_spark_corr) = colnames(mb_spark_corr) = gene_order
  
  write.csv(mb_spark_corr, file.path(dir.outpath, "giotto_results_sim.csv"))
}
