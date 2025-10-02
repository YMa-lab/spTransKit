import argparse
import os
import scanpy as sc
import pandas as pd

def get_r_hvg(norm_counts_dir : str, hvg_dir : str, num_hvg : int):
    files = ["sctransform", "dino", "spanorm"]
    for file in files:
        try:
            matrix = sc.read_h5ad(norm_counts_dir + "/" + file + ".h5ad").to_df()
            hvg = matrix.var(0).sort_values(ascending = False).index.to_list()[0:num_hvg]
            pd.Series(hvg, name = "HVG").to_csv(hvg_dir + "/" + file + ".csv")
        except(FileNotFoundError):
            continue

def merge_runtime_and_memory(results_dir : str):
    runtime_p = pd.read_csv(results_dir + "/runtime_p.csv", header = 0, index_col = 0).iloc[:, 0]
    runtime_r = pd.read_csv(results_dir + "/runtime_r.csv", header = 0, index_col = 0).iloc[:, 0]    
    pd.Series(pd.concat([runtime_p, runtime_r]), name = "Runtime (s)").to_csv(results_dir + "/runtime.csv")
    os.remove(results_dir + "/runtime_p.csv")
    os.remove(results_dir + "/runtime_r.csv")

    memory_p = pd.read_csv(results_dir + "/memory_p.csv", header = 0, index_col = 0).iloc[:, 0]
    memory_r = pd.read_csv(results_dir + "/memory_r.csv", header = 0, index_col = 0).iloc[:, 0]    
    pd.Series(pd.concat([memory_p, memory_r]), name = "Memory (MB)").to_csv(results_dir + "/memory.csv")
    os.remove(results_dir + "/memory_p.csv")
    os.remove(results_dir + "/memory_r.csv")

if __name__ == "__main__":
    # Parse user input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_id", type = str, help = "sample identifier")
    parser.add_argument("--technology", type = str, help = "type of ST technology")
    parser.add_argument("--num_hvg", type = int, help = "number of highly variable genes to be selected")
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

    results_dir = dir + "1.Evaluations/" + params.technology + "/" + params.sample_id + "/july_results"
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    get_r_hvg(norm_counts_dir, hvg_dir, params.num_hvg)
    merge_runtime_and_memory(results_dir)

    print("Transformations completed for " + params.technology + " " + params.sample_id)