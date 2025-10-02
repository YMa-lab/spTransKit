library(argparse)
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

# Get data
coord <- as.matrix(read.csv(paste(c(wd, args[2], args[1], "coordinates.csv"), collapse = "/"), row.names = 1))
x <- t(as.matrix(read_h5ad(paste(c(wd, args[2], args[1], "norm_counts/raw.h5ad"), collapse = "/"))))
row.names(coord) <- colnames(x)
genes <- read.csv("/oscar/data/yma16/Project/spTransform/code/Experiments/data/marker_genes.csv")

# Create subpanel
gene_list <- as.character(subset(genes, Technology == args[2] & Tissue == args[1])[, c("Gene.1", "Gene.2", "Gene.3", "Gene.4", "Gene.5", "Gene.6")])
gene_list <- gene_list[nzchar(gene_list)]
subpanel <- sample(row.names(x), floor(nrow(x) / 4))
for (j in 1:length(gene_list)) {
    if (!(gene_list[j] %in% subpanel)) {
        subpanel <- c(subpanel, gene_list[j])
        subpanel <- subpanel[-1]
    }
}
x <- x[subpanel, ]
x <- x[rowSums(x) != 0, ]
x <- x[, colSums(x) != 0]

coord <- coord[colnames(x), ]

# Create subcounts dir
save_dir <- paste(c(wd, args[2], args[1], "r_subcounts"), collapse = "/")
if (!dir.exists(save_dir)) {
    dir.create(save_dir)
}

# SCTransform
success <- FALSE
counter <- 1
while (!success & counter <= 5) {
    
    trial <- try(vst(x, min_cells = 1, verbosity = 0)[["y"]], silent = TRUE)

    if (!inherits(trial, "try-error")) {
        seurat <- trial
        success <- TRUE
    } else {
        counter <- counter + 1
        Sys.sleep(3)
    }
}

if (success) {
    write_h5ad(anndata = AnnData(t(seurat)), filename = paste(c(save_dir, "sctransform.h5ad"), collapse = "/"), compression = "gzip")
}
    
# Dino
success <- FALSE
counter <- 1
while (!success & counter <= 5) {
    
    trial <- try(Dino(x), silent = TRUE)

    if (!inherits(trial, "try-error")) {
        dino <- trial
        success <- TRUE
    } else {
        counter <- counter + 1
        Sys.sleep(3)
    }
}

if (success) {
    write_h5ad(anndata = AnnData(t(as.matrix(dino))), filename = paste(c(save_dir, "dino.h5ad"), collapse = "/"), compression = "gzip")
}

# SpaNorm
spe <- SpatialExperiment(assays = list(counts = x), spatialCoords = coord)

success <- FALSE
counter <- 1
while (!success & counter <= 5) {
    
    trial <- try(SpaNorm(spe, verbose = FALSE), silent = TRUE)

    if (!inherits(trial, "try-error")) {
        spanorm <- trial
        success <- TRUE
    } else {
        counter <- counter + 1
        Sys.sleep(3)
    }
}

if (success) {
    write_h5ad(anndata = AnnData(t(as.matrix(spanorm@assays@data@listData[["logcounts"]]))), filename = paste(c(save_dir, "spanorm.h5ad"), collapse = "/"), compression = "gzip")
}
