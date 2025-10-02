library(ggplot2)
library(heatmaply)
library(dplyr)
library(htmlwidgets)
library(tidyr)

technology <- 'seqFISH+'
sample.name <- 'Mouse_SVZ'
sample.name2 <- 'Mouse_SVZ'

print(paste("Sample Name started:", sample.name))

dir.output <- file.path('/oscar/data/yma16/Project/spTransform/2.Downstream/results', technology, technology, sample.name, 'Giotto_SCG_july')
filepath_norm <- file.path('/oscar/data/yma16/Project/spTransform/0.NormalizedCounts/', technology, sample.name2,'norm_counts')
coord_file <- file.path('/oscar/data/yma16/Project/spTransform/0.NormalizedCounts', technology, sample.name2, 'coordinates.csv')

coord <- read.csv(coord_file, row.names = 1)
norm_list <- list.files(filepath_norm, full.names = TRUE)
norm_list <- norm_list[file.info(norm_list)$isdir == FALSE]
norm_list <- norm_list[grepl("\\.h5ad$", norm_list)]
norm_list <- basename(norm_list)
norm_list <- sub("\\.h5ad$", "", norm_list)
norm_list <- c(norm_list, "original_giotto")
norm_list <- norm_list[! norm_list %in% c("scanpy_seurat", "scanpy_weinreb", "deseq2_log1p")]

name_dict <- c(
  "raw" = "Raw Counts w/o Trans.",
  "raw_size" = "y/s",
  "cpm" = "CPM",
  "shifted_log" = "log(y/s + 1)",
  "cpm_shifted_log" = "log(CPM + 1)",
  "shifted_log_size" = "log(y/s + 1)/u",
  "acosh" = "acosh(2αy/s + 1)",
  "pseudo_shifted_log" = "log(y/s + 1/4α)",
  "analytic_pearson_residual_clip" = "Analytic Pearson (clip)",
  "analytic_pearson_residual_noclip" = "Analytic Pearson (no clip)",
  "scanpy_zheng" = "scanpy Zheng",
  "scanpy_pearson_residual" = "scanpy Pearson Residual",
  "deseq2_log" = "DESeq2 (log)",
  "deseq2_log1p" = "DESeq2 (log1p)",
  "dino" = "Dino",
  "tmm" = "TMM",
  "normalisr" = "Normalisr",
  "psinorm" = "PsiNorm",
  "original_giotto" = "Original Giotto", 
  "sctransform" = "SCTransform ", 
  "spanorm" = "SpaNorm"
)

all_sig_genes <- character()

for (normalization in norm_list) {
  print(paste("Loading data for:", normalization))
  dir.outpath <- file.path(dir.output, normalization)
  
  km_spatialgenes <- read.csv(file.path(dir.outpath, "giotto_scg_kmeans_binspect.csv"))
  print(sum(km_spatialgenes$adj.p.value < 0.05, na.rm = TRUE))
  print(sum(km_spatialgenes$estimate > 2, na.rm = TRUE))

  keep <- with(km_spatialgenes, adj.p.value < 0.05 & estimate > 2)
  genes <- unique(km_spatialgenes$gene[keep])
  genes <- genes[!is.na(genes) & genes != ""]           # clean NAs/empties
  all_sig_genes <- union(all_sig_genes, genes)          # set-union accumulate
}
sig_gene_set <- sort(all_sig_genes)  # final “set” (unique vector)
cat("Total unique genes:", length(sig_gene_set), "\n")

write.csv(data.frame(gene = sig_gene_set),
          "/oscar/data/yma16/Project/spTransform/2.Downstream/results/summary/SCG/all_spatial_genes_new.csv",
          row.names = FALSE)

km_score_df <- data.frame()
for (normalization in norm_list) {
  print(paste("Loading data for:", normalization))
  dir.outpath <- file.path(dir.output, normalization)
  
  km_spatialgenes <- read.csv(file.path(dir.outpath, "giotto_scg_kmeans_binspect.csv"))[, c("genes", "estimate")]
  km_spatialgenes_subset <- km_spatialgenes[km_spatialgenes$gene %in% sig_gene_set,]
  print(dim(km_spatialgenes_subset))

  colnames(km_spatialgenes_subset)[2] <-normalization
  
  if (nrow(km_score_df) == 0) {
    km_score_df <- km_spatialgenes_subset
  } else {
    km_score_df <- merge(km_score_df, km_spatialgenes_subset, by = "genes", all = TRUE)
  }
}

rownames(km_score_df) <- km_score_df[[1]]   # set first column as row names
km_score_df <- km_score_df[, -1, drop = FALSE]
km_score_df <- km_score_df[, setdiff(names(km_score_df),
                                     c("raw", "original_giotto")),
                           drop = FALSE]
data_frame <- km_score_df
X <- as.matrix(data_frame)
n <- nrow(X); p <- ncol(X)
R01 <- matrix(NA_real_, nrow = n, ncol = p,
              dimnames = list(rownames(X), colnames(X)))
for (i in seq_len(n)) {
  v <- as.numeric(X[i, ])                      # vector for gene i (length = p)
  r <- rank(v, ties.method = "average")        # ranks within the gene
  R01[i, ] <- if (p > 1) (r - 1) / (p - 1) else 0
}
rank_mean <- colMeans(R01, na.rm = TRUE)
rank_norm <- (rank(rank_mean, ties.method = "average") - 1) / (length(rank_mean) - 1)
out_row <- as.data.frame(t(rank_norm))   # 1 row, 20 columns (keeps names as headers)
write.csv(out_row, "../../data/SCG_rank_new_v3.csv", row.names = FALSE)
rank_mean_out <- as.data.frame(t(R01))
write.csv(rank_mean_out, "/oscar/data/yma16/Project/spTransform/2.Downstream/results/summary/SCG/SCG_rank_new_all_v3.csv", row.names = TRUE)

