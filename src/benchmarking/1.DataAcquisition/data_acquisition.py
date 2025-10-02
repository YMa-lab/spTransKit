import argparse
import os
import scanpy as sc

import da_helpers as h

def acquire_data(params : argparse.Namespace, coordinates_dir : str, norm_counts_dir : str):
    # Handle gene count and coordinate data
    coordinates = h.handle_coordinates(params.coordinates, params.technology)
    counts = h.handle_counts(params.counts, coordinates, params.technology)
    # Check to make sure matrices have same dimensions
    coordinates = coordinates.loc[counts.index.to_list(), :]
    h.check_dimensions(counts, coordinates)
    # Save counts and coordinates
    coordinates.to_csv(coordinates_dir + "/coordinates.csv")
    sc.AnnData(counts).write_h5ad(norm_counts_dir + "/raw.h5ad", compression = "gzip")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_id", type = str, help = "sample identifier")
    parser.add_argument("--technology", type = str, help = "type of ST technology")
    parser.add_argument("--counts", type = str, help = "file with gene counts")
    parser.add_argument("--coordinates", type = str, help = "file with spot coordinates")
    params = parser.parse_args()

    # Create directories
    dir = "/oscar/data/yma16/Project/spTransform/"

    coordinates_dir = dir + "0.NormalizedCounts/" + params.technology + "/" + params.sample_id
    if not os.path.exists(coordinates_dir):
        os.makedirs(coordinates_dir)

    norm_counts_dir = coordinates_dir + "/norm_counts"
    if not os.path.exists(norm_counts_dir):
        os.makedirs(norm_counts_dir)

    # 1. Data Acquisition and Preprocessing
    acquire_data(params, coordinates_dir, norm_counts_dir)

    print("Data acquired for " + params.technology + " " + params.sample_id)