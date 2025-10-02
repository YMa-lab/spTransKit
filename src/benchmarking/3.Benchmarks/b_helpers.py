import pandas as pd
import scanpy as sc
import numpy as np
import scipy.sparse

def handle_marker_genes(technology : str, sample : str) -> dict:
    x = pd.read_csv("Experiments/data/marker_genes.csv", header = 0)

    row = x[x["Tissue"] == sample]
    row = pd.Series(row.values.tolist()[0], index = row.columns)
    is_na = row.isna()
    
    genes = {}
    genes[row.loc["Gene 1"]] = row.loc["Distribution"]
    genes[row.loc["Gene 2"]] = row.loc["Distribution.1"]
    genes[row.loc["Gene 3"]] = row.loc["Distribution.2"]
    if is_na.loc["Gene 4"] == False:
        genes[row.loc["Gene 4"]] = row.loc["Distribution.3"]
    if is_na.loc["Gene 5"] == False:
        genes[row.loc["Gene 5"]] = row.loc["Distribution.4"]
    if is_na.loc["Gene 6"] == False:
        genes[row.loc["Gene 6"]] = row.loc["Distribution.5"]
    
    return genes

def handle_marker_pathways(technology : str, sample : str) -> list:
    x = pd.read_csv("Experiments/data/marker_pathways.csv", header = 0)
    row = x[(x["Technology"] == technology) & (x["Tissue"] == sample)].values.tolist()[0][3:]
    
    return [x for x in row if str(x) != "nan"]

def read_counts(files : list, norm_counts_dir : str) -> list[pd.DataFrame]:
    transformations = []
    for file in files:
        try:
            transformations.append(sc.read_h5ad(norm_counts_dir + "/" + file + ".h5ad").to_df())
        except(FileNotFoundError):
            transformations.append(None)
    return transformations

def read_hvg(files : list, num_hvg : int, hvg_dir : str) -> list[list]:
    hvgs = []
    for file in files:
        try:
            hvgs.append(pd.read_csv(hvg_dir + "/" + file + ".csv", header = 0, index_col = 0).loc[:, "HVG"].tolist()[0:num_hvg])
        except(FileNotFoundError):
            hvgs.append(None)
    return hvgs

def weight_matrix(coordinates : pd.DataFrame) -> pd.DataFrame:
    length = len(coordinates.index)
    w = scipy.sparse.dok_matrix((length, length)).astype(np.int64)
    for i in range(0, len(coordinates.index)):
        ix = coordinates.iloc[i, 0]
        iy = coordinates.iloc[i, 1]
        
        xdist = (coordinates.loc[:, "x"] - ix) ** 2
        ydist = (coordinates.loc[:, "y"] - iy) ** 2
        dist = pd.DataFrame(np.exp(-((xdist + ydist) ** 0.5)), index = coordinates.index)
        dist["pos"] = list(range(0, len(coordinates.index)))

        dist = dist.loc[dist[0] < 1.0, :]
        dist = dist.loc[dist[0] >= dist.loc[:, 0].mean(), :]
        for j in dist.index:
            w[i, dist.loc[j, "pos"]] = np.ceil(dist.loc[j, 0])

    return w

def morans_i(weight_matrix : scipy.sparse.dok_matrix, gene_dist : pd.Series) -> float:
    n = gene_dist.size
    mean = gene_dist.mean()
    i_dist = gene_dist - mean
    i_dist_squared = i_dist ** 2
    j_dist = gene_dist - mean
    numerator = weight_matrix.copy()
    numerator = numerator.multiply(i_dist)
    numerator = numerator.T.multiply(j_dist)
    weight_sum = weight_matrix.sum()
    if weight_sum * i_dist_squared.sum() == 0.0:
        moran = 0.0
    else:
        moran = (n * numerator.sum(0).sum()) / (weight_sum * i_dist_squared.sum())
    return moran

def rank(results : list[pd.Series], names : list) -> list[pd.Series]:
    new_results = []

    # HVG Similarity
    x = results[0].loc[~results[0].index.isin(["y"])]
    x.sort_values(ascending = False, inplace = True)
    n = len(x.dropna())
    x[0:n] = list(range(1, n + 1))
    x["y"] = 0
    x = x.reindex(names)
    new_results.append(x)

    # HVG Consistency
    x = results[1].loc[~results[1].index.isin(["y"])]
    x.sort_values(ascending = False, inplace = True)
    n = len(x.dropna())
    x[0:n] = list(range(1, n + 1))
    x["y"] = 0
    x = x.reindex(names)
    new_results.append(x)

    # Gene Gene Correlation
    x = results[2].loc[~results[2].index.isin(["y"])]
    x.sort_values(ascending = True, inplace = True)
    n = len(x.dropna())
    x[0:n] = list(range(1, n + 1))
    x["y"] = 0
    x = x.reindex(names)
    new_results.append(x)

    # k NN Similarity
    x = results[3].loc[~results[3].index.isin(["y"])]
    x.sort_values(ascending = False, inplace = True)
    n = len(x.dropna())
    x[0:n] = list(range(1, n + 1))
    x["y"] = 0
    x = x.reindex(names)
    new_results.append(x)

    # Marker Gene Distribution
    x = results[4].loc[~results[4].index.isin(["y"])]
    x.sort_values(ascending = True, inplace = True)
    n = len(x.dropna())
    x[0:n] = list(range(1, n + 1))
    x["y"] = 0
    x = x.reindex(names)
    new_results.append(x)

    # Pathway Enrichment
    x = results[5].loc[~results[5].index.isin(["y"])]
    x.sort_values(ascending = False, inplace = True)
    n = len(x.dropna())
    x[0:n] = list(range(1, n + 1))
    x["y"] = 0
    x = x.reindex(names)
    new_results.append(x)

    # Runtime
    x = results[6].loc[~results[6].index.isin(["y"])]
    x.sort_values(ascending = True, inplace = True)
    n = len(x.dropna())
    x[0:n] = list(range(1, n + 1))
    x["y"] = 0
    x = x.reindex(names)
    new_results.append(x)

    # Memory
    x = results[7].loc[~results[7].index.isin(["y"])]
    x.sort_values(ascending = True, inplace = True)
    n = len(x.dropna())
    x[0:n] = list(range(1, n + 1))
    x["y"] = 0
    x = x.reindex(names)
    new_results.append(x)

    # Size Factor Correspondence
    x = results[8].loc[~results[8].index.isin(["y"])]
    x.sort_values(ascending = False, inplace = True)
    n = len(x.dropna())
    x[0:n] = list(range(1, n + 1))
    x["y"] = 0
    x = x.reindex(names)
    new_results.append(x)

    # Count Correspondence
    x = results[9].loc[~results[9].index.isin(["y"])]
    x.sort_values(ascending = False, inplace = True)
    n = len(x.dropna())
    x[0:n] = list(range(1, n + 1))
    x["y"] = 0
    x = x.reindex(names)
    new_results.append(x)

    return new_results
