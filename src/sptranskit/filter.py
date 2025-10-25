import scanpy as sc
import pandas as pd


def filter_counts(data : sc.AnnData,
                  min_counts_per_gene : int = 1,
                  min_counts_per_loc : int = 1,
                  inplace : bool = False):
    """ This function filters the gene count matrix to remove low quality spatial locations and genes.
    
    Parameters
    ----------
    data
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
    min_counts_per_gene
        int, the minimum number of counts a gene must have across all spatial locations, default = 1.
    min_counts_per_loc
        int, the minimum number of counts a spatial location must have across all genes, default = 1.
    inplace
        bool, whether to filter the original AnnData object, default = False.

    Returns
    ----------
    copy
        sc.AnnData, copy of the filtered gene count matrix if inplace = False.

    """

    if data.obsm["spatial"] is None:
        raise Exception("No spatial information included in AnnData object.")
    
    counts = data.to_df()
    counts = counts.loc[:, counts.sum(0) >= min_counts_per_gene]
    counts = counts.loc[counts.sum(1) >= min_counts_per_loc, :]
    counts = counts.loc[:, ~counts.columns.duplicated()]
    counts = counts.loc[~counts.index.duplicated(), :]
    
    if inplace:
        data.var_names_make_unique()
        data = data[counts.index.to_list(), counts.columns.to_list()]
    else:
        copy = data.copy()
        copy.var_names_make_unique()
        copy = copy[counts.index.to_list(), counts.columns.to_list()]
        return copy
