library(ggplot2)
library(heatmaply)
library(dplyr)
library(htmlwidgets)
library(tidyr)

# font_add("Arial", regular = "/oscar/data/yma16/Project/spTransform/2.Downstream/Arial Unicode.ttf")  # use actual path on cluster
# showtext_auto()

technology <- 'seqFISH+'
sample.name <- 'Mouse_SVZ'
sample.name2 <- 'Mouse_SVZ'

print(paste("Sample Name started:", sample.name))

dir.output <- file.path('/oscar/data/yma16/Project/spTransform/2.Downstream/results', technology, technology, sample.name, 'Giotto_SCG_july_set')
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


############### plot correlation and RMSE plot ######################
# Initialize list to store flattened matrices for all normalizations
flattened_matrices <- list()
norm_methods <- norm_list

# Step 1: Load each correlation matrix, flatten it, and store it
ref_genes <- NULL
if (length(norm_list) > 0) {
  ref_genes <- rownames(read.csv(file.path(dir.output, norm_list[1], "giotto_results_sim.csv"), row.names = 1))
}

for (normalization in norm_list) {
  dir.outpath <- file.path(dir.output, normalization)
  mat <- read.csv(file.path(dir.outpath, "giotto_results_sim.csv"), row.names = 1, check.names = FALSE)

  # mat <- mat %>% select_if(is.numeric)
  # mat[] <- lapply(mat, function(x) if(is.factor(x)) as.numeric(as.character(x)) else as.numeric(x))
  mat <- as.matrix(mat)

  # colnames(mat) <- gsub("^X", "", colnames(mat))
  mat <- mat[ref_genes, ref_genes, drop = FALSE]

  # Flatten the matrix by extracting the lower triangle (excluding diagonal) to avoid redundancy
  flat_vec <- as.vector(mat[lower.tri(mat)])
  flattened_matrices[[normalization]] <- flat_vec
}

# Step 2: Compute pairwise Pearson correlation and RMSE for each normalization method
cor_results <- data.frame()
rmse_results <- data.frame()

# Only include normalization methods with valid flattened matrices
valid_norm_methods <- names(flattened_matrices)

for (i in 1:length(valid_norm_methods)) {
  for (j in (i+1):length(valid_norm_methods)) {
    norm1 <- valid_norm_methods[i]
    norm2 <- valid_norm_methods[j]
    
    # Flattened matrices
    flat1 <- flattened_matrices[[norm1]]
    flat2 <- flattened_matrices[[norm2]]
    
    # Ensure lengths match and remove NA values
    if (length(flat1) == length(flat2)) {
      # Remove any NA values in flat1 and flat2
      na_mask <- !is.na(flat1) & !is.na(flat2)
      flat1 <- flat1[na_mask]
      flat2 <- flat2[na_mask]
      
      # Check length again after removing NA values
      if (length(flat1) == length(flat2) && length(flat1) > 0) {
        # Calculate Pearson correlation
        pearson_corr <- cor(flat1, flat2, method = "pearson")
        
        # Calculate RMSE
        rmse <- sqrt(mean((flat1 - flat2)^2))
        
        # Store results
        cor_results <- rbind(cor_results, data.frame(Norm1 = norm1, Norm2 = norm2, Correlation = pearson_corr))
        rmse_results <- rbind(rmse_results, data.frame(Norm1 = norm1, Norm2 = norm2, RMSE = rmse))
      } else {
        warning(paste("Length mismatch or empty vectors after NA removal for", norm1, "and", norm2))
      }
    } else {
      warning(paste("Initial length mismatch for", norm1, "and", norm2))
    }
  }
}

# Step 3: Prepare data for plotting
# Reshape correlation and RMSE results to long format for ggplot
cor_long <- cor_results %>%
  pivot_longer(cols = c("Norm1", "Norm2"), names_to = "Metric", values_to = "Normalization") %>%
  select(-Metric) %>%
  rename(Value = Correlation)

rmse_long <- rmse_results %>%
  pivot_longer(cols = c("Norm1", "Norm2"), names_to = "Metric", values_to = "Normalization") %>%
  select(-Metric) %>%
  rename(Value = RMSE)

write.csv(cor_long, "/oscar/data/yma16/Project/spTransform/2.Downstream/results/summary/SCE/cor_all_scores_v2.csv", row.names = FALSE)
write.csv(rmse_long, "/oscar/data/yma16/Project/spTransform/2.Downstream/results/summary/SCE/rmse_all_scores_v2.csv", row.names = FALSE)

generate_normalized_rank_list <- function(rmse_long, cor_long) {
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
  
  # Median RMSE ranking
  rmse_ranks <- rmse_long %>%
    group_by(Normalization) %>%
    summarise(MedianRMSE = median(Value, na.rm = TRUE)) %>%
    filter(!Normalization %in% c("Raw Counts w/o Trans.", "Original Giotto")) %>%
    mutate(Rank_RMSE = rank(-MedianRMSE, ties.method = "average"),
           NormalizedRank_RMSE = (Rank_RMSE - 1) / (max(Rank_RMSE) - 1))
  
  # Median correlation ranking
  cor_ranks <- cor_long %>%
    group_by(Normalization) %>%
    summarise(MedianCorrelation = median(Value, na.rm = TRUE)) %>%
    filter(!Normalization %in% c("Raw Counts w/o Trans.", "Original Giotto")) %>%
    mutate(Rank_Correlation = rank(MedianCorrelation, ties.method = "average"),
           NormalizedRank_Correlation = (Rank_Correlation - 1) / (max(Rank_Correlation) - 1))
  
  # Combine the RMSE and correlation rankings
  final_ranks <- left_join(rmse_ranks[, c("Normalization", "NormalizedRank_RMSE")],
                           cor_ranks[, c("Normalization", "NormalizedRank_Correlation")],
                           by = "Normalization")
  
  # Compute the average rank
  final_ranks$AverageRank <- rowMeans(final_ranks[, c("NormalizedRank_RMSE", "NormalizedRank_Correlation")])

  
  # Save the final rank list to a CSV file
  output_file <- file.path("/oscar/data/yma16/Project/spTransform/2.Downstream/results/summary/SCE/SCE_all_tests_normalized_ranks_v2.csv")
  write.csv(final_ranks, file = output_file, row.names = FALSE)
  write.csv(final_ranks[, c("Normalization", "AverageRank")], file = file.path("/oscar/data/yma16/Project/spTransform/2.Downstream/results/summary/SCE/SCE_total_rank_v2.csv"), row.names = FALSE)
  write.csv(final_ranks[, c("Normalization", "AverageRank")], file = file.path("/oscar/data/yma16/Project/spTransform/code/Manuscript/Version5/data/SCE_total_rank_v2.csv"), row.names = FALSE)
  
  return(final_ranks)
}

final_ranks <- generate_normalized_rank_list(rmse_long, cor_long) 
