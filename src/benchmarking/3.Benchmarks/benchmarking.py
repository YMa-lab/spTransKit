import argparse
import os
import pandas as pd

import b_helpers as h
import bio_benchmarks as b
import nonbio_benchmarks as n

def benchmark(params : argparse.Namespace, names : list, files : list, norm_counts_dir : str, coordinates_dir : str, hvg_dir : str, results_dir : str):
    # Handle marker genes and marker pathways
    genes = h.handle_marker_genes(params.technology, params.sample_id)
    pathways = h.handle_marker_pathways(params.technology, params.sample_id)

    # Read files
    trans = h.read_counts(files, norm_counts_dir)
    coordinates = pd.read_csv(coordinates_dir + "/coordinates.csv", header = 0, index_col = 0)
    hvgs = h.read_hvg(files, params.num_hvg, hvg_dir)

    # Run biological benchmarks
    hvg_sim, hvg_cons = b.hvg_similarity(hvgs, names, results_dir)
    gene_corr = b.gene_correlation(trans, hvgs, names, results_dir)
    knn_sim = b.knn_similarity(trans, hvgs, files, names, params.num_neighbors, knn_dir, results_dir)
    gene_dist = b.gene_distributions(coordinates, genes, trans, names, results_dir)
    pwy_enr = b.pathway_enrichment(pathways, hvgs, names, params.species, results_dir)

    # Run non-biological benchmarks
    runtime = n.get_runtime(results_dir, names)
    memory = n.get_memory(results_dir, names)
    subpanel = n.create_subpanel(trans[0], genes)
    t_two = n.transform_again(trans[0], names, params.num_hvg, subpanel, coordinates_dir)
    sf_corr = n.sf_correspondence(trans, t_two, names, results_dir)
    count_corr = n.count_correspondence(trans, t_two, names, genes, results_dir)

    # Save total results
    rank = pd.concat(h.rank([hvg_sim, hvg_cons, gene_corr, knn_sim, gene_dist, pwy_enr, runtime, memory, sf_corr, count_corr], names), axis = 1)
    rank.columns = ["HVG Similarity", "HVG Consistency", "Gene Gene Correlation", "k NN Similarity", "Marker Gene Distribution", "Pathway Enrichment", "Runtime", "Memory", "Size Factor Correspondence", "Count Correspondence"]
    rank.to_csv(results_dir + "/rank_matrix.csv")

if __name__ == "__main__":
    # Parse user input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_id", type = str, help = "sample identifier")
    parser.add_argument("--technology", type = str, help = "type of ST technology")
    parser.add_argument("--species", type = str, help = "tissue species")
    parser.add_argument("--num_hvg", type = int, help = "number of highly variable genes to be selected")
    parser.add_argument("--num_neighbors", type = int, help = "number of nearest neighbors to be selected")
    params = parser.parse_args()

    # Create directories
    dir = "/oscar/data/yma16/Project/spTransform/"

    coordinates_dir = dir + "0.NormalizedCounts/" + params.technology + "/" + params.sample_id
    if not os.path.exists(coordinates_dir):
        os.makedirs(coordinates_dir)

    norm_counts_dir = coordinates_dir + "/norm_counts"
    if not os.path.exists(norm_counts_dir):
        os.makedirs(norm_counts_dir)

    hvg_dir = dir + "1.Evaluations/" + params.technology + "/" + params.sample_id + "/july_hvg"
    if not os.path.exists(hvg_dir):
        os.makedirs(hvg_dir)

    knn_dir = dir + "1.Evaluations/" + params.technology + "/" + params.sample_id + "/knn"
    if not os.path.exists(knn_dir):
        os.makedirs(knn_dir)

    results_dir = dir + "1.Evaluations/" + params.technology + "/" + params.sample_id + "/july_results"
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    # Setup data
    names = ["y", "y/s", "CPM", "log(y/s + 1)", "acosh(2αy/s + 1)", "log(y/s + 1/(4α))", "log(CPM + 1)", "log(y/s + 1)/u", "Analytic Pearson (no clip)", "Analytic Pearson (clip)", "scanpy Zheng", "scanpy Pearson Residual", "DESeq2 (log)", "DESeq2 (log1p)", "TMM", "Normalisr", "PsiNorm", "SCTransform", "Dino", "SpaNorm"]
    files = ["raw", "raw_size", "cpm", "shifted_log", "acosh", "pseudo_shifted_log", "cpm_shifted_log", "shifted_log_size", "analytic_pearson_residual_noclip", "analytic_pearson_residual_clip", "scanpy_zheng", "scanpy_pearson_residual", "deseq2_log", "deseq2_log1p", "tmm", "normalisr", "psinorm", "sctransform", "dino", "spanorm"]
    
    # 3. Benchmarking Analysis
    benchmark(params, names, files, norm_counts_dir, coordinates_dir, hvg_dir, results_dir)

    print("Benchmarks completed for " + params.technology + " " + params.sample_id)
    