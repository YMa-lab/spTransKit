import scanpy as sc
import pandas as pd


if __name__ == "__main__":
    matrix = sc.read_10x_h5("src/sptranskit/data/DLPFC_151673/unfiltered_raw.h5")
    coordinates = pd.read_csv("src/sptranskit/data/DLPFC_151673/unfiltered_coord.txt", header = None, index_col = 0)

    coordinates = coordinates.loc[coordinates[1] != 0, 2:3]
    coordinates.columns = ["x", "y"]

    data = sc.AnnData(matrix.to_df().loc[coordinates.index, :])
    data.obsm["spatial"] = coordinates

    counts = data.to_df()
    counts = counts.loc[:, counts.sum(0) >= 1]
    counts = counts.loc[counts.sum(1) >= 1, :]
    counts = counts.loc[:, ~counts.columns.duplicated()]
    counts = counts.loc[~counts.index.duplicated(), :]

    data.var_names_make_unique()
    data = data[counts.index.to_list(), counts.columns.to_list()]
    data.raw = data.copy()
    print(data.raw)
