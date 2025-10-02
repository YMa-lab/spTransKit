library(argparse)
library(pryr)
library(anndata)
library(sctransform)
library(Dino)
library(SpaNorm)
library(SpatialExperiment)

parser <- ArgumentParser()
parser$add_argument("--sample_id", type = "character")
parser$add_argument("--technology", type = "character")
args <- parser$parse_args()

wd <- "/oscar/data/yma16/Project/spTransform/0.NormalizedCounts"

# Runtime and memory
r <- data.frame(matrix(nrow = 3, ncol = 1))
row.names(r) <- c("SCTransform", "Dino", "SpaNorm")
colnames(r) <- c("Runtime")
m <- data.frame(matrix(nrow = 3, ncol = 1))
row.names(m) <- c("SCTransform", "Dino", "SpaNorm")
colnames(m) <- c("Memory")
    
# Get data
coord <- as.matrix(read.csv(paste(c(wd, args[2], args[1], "coordinates.csv"), collapse = "/"), row.names = 1))
x <- t(as.matrix(read_h5ad(paste(c(wd, args[2], args[1], "norm_counts/raw.h5ad"), collapse = "/"))))

# SCTransform
success <- FALSE
counter <- 1
while (!success & counter <= 5) {
    
    start_run <- Sys.time()
    start_mem <- mem_used()
    trial <- try(vst(x, min_cells = 1, verbosity = 0)[["y"]], silent = TRUE)
    end_mem <- mem_used()
    end_run <- Sys.time()

    if (!inherits(trial, "try-error")) {
        seurat <- trial
        success <- TRUE
    } else {
        counter <- counter + 1
        Sys.sleep(3)
    }
}
    
if (success) {
    r["SCTransform", "Runtime"] <- end_run - start_run
    m["SCTransform", "Memory"] <- (end_mem - start_mem) / 1000000
    write_h5ad(anndata = AnnData(t(seurat)), filename = paste(c(wd, args[2], args[1], "norm_counts/sctransform.h5ad"), collapse = "/"), compression = "gzip")
} else {
    r["SCTransform", "Runtime"] <- NA
    m["SCTransform", "Memory"] <- NA
}

# Dino
success <- FALSE
counter <- 1
while (!success & counter <= 5) {
    
    start_run <- Sys.time()
    start_mem <- mem_used()
    trial <- try(Dino(x), silent = TRUE)
    end_mem <- mem_used()
    end_run <- Sys.time()

    if (!inherits(trial, "try-error")) {
        dino <- trial
        success <- TRUE
    } else {
        counter <- counter + 1
        Sys.sleep(3)
    }
}

if (success) {
    r["Dino", "Runtime"] <- end_run - start_run
    m["Dino", "Memory"] <- (end_mem - start_mem) / 1000000
    write_h5ad(anndata = AnnData(t(as.matrix(dino))), filename = paste(c(wd, args[2], args[1], "norm_counts/dino.h5ad"), collapse = "/"), compression = "gzip")
} else {
    r["Dino", "Runtime"] <- NA
    m["Dino", "Memory"] <- NA
}

# SpaNorm
spe <- SpatialExperiment(assays = list(counts = x), spatialCoords = coord)

success <- FALSE
counter <- 1
while (!success & counter <= 5) {
    
    start_run <- Sys.time()
    start_mem <- mem_used()
    trial <- try(SpaNorm(spe, verbose = FALSE), silent = TRUE)
    end_mem <- mem_used()
    end_run <- Sys.time()

    if (!inherits(trial, "try-error")) {
        spanorm <- trial
        success <- TRUE
    } else {
        counter <- counter + 1
        Sys.sleep(3)
    }
}

if (success) {
    r["SpaNorm", "Runtime"] <- end_run - start_run
    m["SpaNorm", "Memory"] <- (end_mem - start_mem) / 1000000
    write_h5ad(anndata = AnnData(t(as.matrix(spanorm@assays@data@listData[["logcounts"]]))), filename = paste(c(wd, args[2], args[1], "norm_counts/spanorm.h5ad"), collapse = "/"), compression = "gzip")
} else {
    r["SpaNorm", "Runtime"] <- NA
    m["SpaNorm", "Memory"] <- NA
}

# Write runtime and memory files
write.csv(r, paste(c(gsub("0.NormalizedCounts", "1.Evaluations", wd), args[2], args[1], "july_results/runtime_r.csv"), collapse = "/"))
write.csv(m, paste(c(gsub("0.NormalizedCounts", "1.Evaluations", wd), args[2], args[1], "july_results/memory_r.csv"), collapse = "/"))
