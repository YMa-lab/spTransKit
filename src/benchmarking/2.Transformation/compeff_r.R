library(pryr)
library(anndata)
library(sctransform)
library(Dino)
library(SpaNorm)
library(SpatialExperiment)

compeff_spot <- function(x, coord, number) {
    # Subset
    sub <- sample(colnames(x), number)
    x <- x[, sub]
    x <- x[rowSums(x) != 0, ]
    x <- x[, colSums(x) != 0]
    coord <- coord[colnames(x), ]
    
    # Runtime and memory
    r <- data.frame(matrix(nrow = 3, ncol = 1))
    row.names(r) <- c("SCTransform", "Dino", "SpaNorm")
    colnames(r) <- c(as.character(number))
    m <- data.frame(matrix(nrow = 3, ncol = 1))
    row.names(m) <- c("SCTransform", "Dino", "SpaNorm")
    colnames(m) <- c(as.character(number))
    
    # SCTransform
    start_run <- Sys.time()
    start_mem <- mem_used()
    vst(x, min_cells = 1, verbosity = 0)[["y"]]
    end_mem <- mem_used()
    end_run <- Sys.time()
    
    r["SCTransform", as.character(number)] <- end_run - start_run
    m["SCTransform", as.character(number)] <- (end_mem - start_mem) / 1000000
    
    # Dino
    start_run <- Sys.time()
    start_mem <- mem_used()
    Dino(x)
    end_mem <- mem_used()
    end_run <- Sys.time()
    
    r["Dino", as.character(number)] <- end_run - start_run
    m["Dino", as.character(number)] <- (end_mem - start_mem) / 1000000
    
    # SpaNorm
    spe <- SpatialExperiment(assays = list(counts = x), spatialCoords = coord)
    
    start_run <- Sys.time()
    start_mem <- mem_used()
    SpaNorm(spe, verbose = FALSE)
    end_mem <- mem_used()
    end_run <- Sys.time()
    
    r["SpaNorm", as.character(number)] <- end_run - start_run
    m["SpaNorm", as.character(number)] <- (end_mem - start_mem) / 1000000
    
    return(list(runtime = r, memory = m))
}

# Computational efficiency for spots
counts <- t(as.matrix(read_h5ad("/oscar/data/yma16/Project/spTransform/0.NormalizedCounts/10x_Xenium_Prime/Human_Prostate/norm_counts/raw.h5ad")))
coordinates <- as.matrix(read.csv("/oscar/data/yma16/Project/spTransform/0.NormalizedCounts/10x_Xenium_Prime/Human_Prostate/coordinates.csv", row.names = 1))

sub_1000 <- compeff_spot(counts, coordinates, 1000)
sub_5000 <- compeff_spot(counts, coordinates, 5000)
sub_10000 <- compeff_spot(counts, coordinates, 10000)
sub_25000 <- compeff_spot(counts, coordinates, 25000)
sub_50000 <- compeff_spot(counts, coordinates, 50000)
sub_100000 <- compeff_spot(counts, coordinates, 100000)
sub_125000 <- compeff_spot(counts, coordinates, 125000)
sub_150000 <- compeff_spot(counts, coordinates, 150000)
sub_175000 <- compeff_spot(counts, coordinates, 175000)
sub_max <- compeff_spot(counts, coordinates, ncol(counts))

spots_r <- cbind(sub_1000[["runtime"]],
                 sub_5000[["runtime"]],
                 sub_10000[["runtime"]],
                 sub_25000[["runtime"]],
                 sub_50000[["runtime"]],
                 sub_100000[["runtime"]],
                 sub_125000[["runtime"]],
                 sub_150000[["runtime"]],
                 sub_175000[["runtime"]],
                 sub_max[["runtime"]]
                 )
write.csv(spots_r, file = "/oscar/data/yma16/Project/spTransform/1.Evaluations/10x_Xenium_Prime/Human_Prostate/july_results/spot_runtime_r.csv")

spots_m <- cbind(sub_1000[["memory"]],
                 sub_5000[["memory"]],
                 sub_10000[["memory"]],
                 sub_25000[["memory"]],
                 sub_50000[["memory"]],
                 sub_100000[["memory"]],
                 sub_125000[["memory"]],
                 sub_150000[["memory"]],
                 sub_175000[["memory"]],
                 sub_max[["memory"]]
                 )
write.csv(spots_m, file = "/oscar/data/yma16/Project/spTransform/1.Evaluations/10x_Xenium_Prime/Human_Prostate/july_results/spot_memory_r.csv")

compeff_gene <- function(x, coord, number) {
    # Subset
    sub <- sample(row.names(x), number)
    x <- x[sub, ]
    x <- x[rowSums(x) != 0, ]
    x <- x[, colSums(x) != 0]
    coord <- coord[colnames(x), ]
    
    # Runtime and memory
    r <- data.frame(matrix(nrow = 3, ncol = 1))
    row.names(r) <- c("SCTransform", "Dino", "SpaNorm")
    colnames(r) <- c(as.character(number))
    m <- data.frame(matrix(nrow = 3, ncol = 1))
    row.names(m) <- c("SCTransform", "Dino", "SpaNorm")
    colnames(m) <- c(as.character(number))
    
    # SCTransform
    start_run <- Sys.time()
    start_mem <- mem_used()
    vst(x, min_cells = 1, verbosity = 0)[["y"]]
    end_mem <- mem_used()
    end_run <- Sys.time()
    
    r["SCTransform", as.character(number)] <- end_run - start_run
    m["SCTransform", as.character(number)] <- (end_mem - start_mem) / 1000000
    
    # Dino
    start_run <- Sys.time()
    start_mem <- mem_used()
    Dino(x)
    end_mem <- mem_used()
    end_run <- Sys.time()
    
    r["Dino", as.character(number)] <- end_run - start_run
    m["Dino", as.character(number)] <- (end_mem - start_mem) / 1000000
    
    # SpaNorm
    spe <- SpatialExperiment(assays = list(counts = x), spatialCoords = coord)
    
    start_run <- Sys.time()
    start_mem <- mem_used()
    SpaNorm(spe, verbose = FALSE)
    end_mem <- mem_used()
    end_run <- Sys.time()
    
    r["SpaNorm", as.character(number)] <- end_run - start_run
    m["SpaNorm", as.character(number)] <- (end_mem - start_mem) / 1000000

    return(list(runtime = r, memory = m))
}

# Computational efficiency for genes
counts <- t(as.matrix(read_h5ad("/oscar/data/yma16/Project/spTransform/0.NormalizedCounts/Open-ST/Human_MLN_3/norm_counts/raw.h5ad")))
coordinates <- as.matrix(read.csv("/oscar/data/yma16/Project/spTransform/0.NormalizedCounts/Open-ST/Human_MLN_3/coordinates.csv", row.names = 1))

sub_1000 <- compeff_gene(counts, coordinates, 1000)
sub_5000 <- compeff_gene(counts, coordinates, 5000)
sub_10000 <- compeff_gene(counts, coordinates, 10000)
sub_15000 <- compeff_gene(counts, coordinates, 15000)
sub_20000 <- compeff_gene(counts, coordinates, 20000)
sub_25000 <- compeff_gene(counts, coordinates, 25000)
sub_30000 <- compeff_gene(counts, coordinates, 30000)
sub_35000 <- compeff_gene(counts, coordinates, 35000)
sub_max <- compeff_gene(counts, coordinates, nrow(counts))

genes_r <- cbind(sub_1000[["runtime"]],
                 sub_5000[["runtime"]],
                 sub_10000[["runtime"]],
                 sub_15000[["runtime"]],
                 sub_20000[["runtime"]],
                 sub_25000[["runtime"]],
                 sub_30000[["runtime"]],
                 sub_35000[["runtime"]],
                 sub_max[["runtime"]]
)
write.csv(genes_r, file = "/oscar/data/yma16/Project/spTransform/1.Evaluations/Open-ST/Human_MLN_3/july_results/gene_runtime_r.csv")

genes_m <- cbind(sub_1000[["memory"]],
                 sub_5000[["memory"]],
                 sub_10000[["memory"]],
                 sub_15000[["memory"]],
                 sub_20000[["memory"]],
                 sub_25000[["memory"]],
                 sub_30000[["memory"]],
                 sub_35000[["memory"]],
                 sub_max[["memory"]]
)
write.csv(genes_m, file = "/oscar/data/yma16/Project/spTransform/1.Evaluations/Open-ST/Human_MLN_3/july_results/gene_memory_r.csv")
