library(reticulate)
use_condaenv("giotto_env", conda = "/users/yzhu194/.local/share/r-miniconda/bin/conda", required = TRUE)
reticulate::py_config()
library(Giotto)
library(Matrix)
library(SingleCellExperiment)
anndata_py <- import("anndata", convert = FALSE)
# library(anndata)

args <- commandArgs(trailingOnly = TRUE)
technology <- 'seqFISH+'
sample.name <- 'Mouse_SVZ'
sample.name2 <- 'Mouse_SVZ'
hvg_type <- 'july_v2'

print(paste("Technology started:", technology))
print(paste("Sample Name started:", sample.name))

dir.output <- file.path('/oscar/data/yma16/Project/spTransform/2.Downstream/results/seqFISH+', technology, sample.name, paste0('Giotto_SCG_',hvg_type))
coord_file <- file.path('/oscar/data/yma16/Project/spTransform/0.NormalizedCounts', technology, sample.name2, 'coordinates.csv')

coord <- read.csv(coord_file, row.names = 1)

filepath_norm <- file.path('/oscar/data/yma16/Project/spTransform/0.NormalizedCounts/', technology, sample.name, 'norm_counts')
norm_list <- list.files(filepath_norm, full.names = TRUE)
norm_list <- norm_list[file.info(norm_list)$isdir == FALSE]
norm_list <- basename(norm_list)
norm_list <- sub("\\.h5ad$", "", norm_list)
norm_list <- c(norm_list, "original_giotto")
norm_list <- norm_list[!grepl("raw\\.csv$", basename(norm_list))]

for (normalization in norm_list) {
  print(paste("Normalization Method started:", normalization))
  dir.outpath <- file.path(dir.output, normalization)
  if(!dir.exists(file.path(dir.outpath))){
    dir.create(file.path(dir.outpath), recursive = TRUE)
  }
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

  # processing
  if (normalization == "original_giotto") {
    
    my_giotto_object <- filterGiotto(gobject = my_giotto_object, 
                                     expression_threshold = 0.5, 
                                     gene_det_in_min_cells = 20, 
                                     min_det_genes_per_cell = 0)
    
    my_giotto_object <- normalizeGiotto(gobject = my_giotto_object)
  }
  else {
    my_giotto_object <- filterGiotto(gobject = my_giotto_object, 
                                     expression_threshold = 0, 
                                     gene_det_in_min_cells = 0, 
                                     min_det_genes_per_cell = 0)
    my_giotto_object@norm_expr = expr_dense
  }
  
  # create network (required for binSpect methods)
  my_giotto_object = createSpatialNetwork(gobject = my_giotto_object, minimum_k = 2)
  
  # identify genes with a spatial coherent expression profile
  km_spatialgenes = binSpect(my_giotto_object, bin_method = 'kmeans')

  spatGenePlot(my_giotto_object, expression_values = 'normalized', 
               genes = km_spatialgenes[1:5]$genes, point_size = 3,
               point_shape = 'border', point_border_stroke = 0.1, cow_n_col = 2)
  
  write.csv(km_spatialgenes, file.path(dir.outpath, "giotto_scg_kmeans_binspect.csv"))

  # binSpect rank method
  rnk_spatialgenes = binSpect(my_giotto_object, bin_method = 'rank')
  
  spatGenePlot(my_giotto_object, expression_values = 'normalized', 
               genes = rnk_spatialgenes[1:5]$genes, point_size = 3,
               point_shape = 'border', point_border_stroke = 0.1, cow_n_col = 2)
  
  write.csv(rnk_spatialgenes, file.path(dir.outpath, "giotto_scg_ranked_binspect.csv"))
  
  # silhouetteRank method
  silh_spatialgenes = silhouetteRank(my_giotto_object)
  
  spatGenePlot(my_giotto_object, expression_values = 'normalized', 
               genes = silh_spatialgenes[1:5]$genes,  point_size = 3,
               point_shape = 'border', point_border_stroke = 0.1, cow_n_col = 2)
  
  write.csv(silh_spatialgenes, file.path(dir.outpath, "giotto_scg_silhouetteRank.csv"))

}

