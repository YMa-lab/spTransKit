import scanpy as sc
import pandas as pd


def filter_counts(matrix : sc.AnnData,
                  coordinates : pd.DataFrame,
                  min_counts_per_gene : int = 1,
                  min_counts_per_loc : int = 1):
    """ This function filters the gene count matrix to remove low quality spatial locations and genes.
    
    Parameters
    ----------
    matrix
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
    coordinates
        pd.DataFrame, N x 2 coordinate matrix. Rows correspond to spatial locations and columns correspond to the x and
        y spatial coordinates.
    min_counts_per_gene
        int, the minimum number of counts a gene must have across all spatial locations, default = 1
    min_counts_per_loc
        int, the minimum number of counts a spatial location must have across all genes, default = 1

    Returns
    ----------
    copy
        sc.AnnData, copy of the filtered gene count matrix if inplace = False
    coord
        pd.DataFrame, copy of the filtered spatial coordinate matrix if inplace = False

    """
    
    counts = matrix.to_df()
    counts = counts.loc[:, counts.sum(0) >= min_counts_per_gene]
    counts = counts.loc[counts.sum(1) >= min_counts_per_loc, :]
    counts = counts.loc[:, ~counts.columns.duplicated()]
    counts = counts.loc[~counts.index.duplicated(), :]
    
    coord = coordinates.loc[counts.index.to_list(), :]
    
    check_dimensions(counts, coord)

    return sc.AnnData(counts), coord


def check_dimensions(counts : pd.DataFrame, coordinates : pd.DataFrame):
    """ This function checks that the filtered gene count matrix and the filtered spatial coordinate matrix share the
        same number of spatial coordinates to ensure data congruence. This function is designed to be called within the
        structure of the filter_counts function.

        Parameters
        ----------
        matrix
            sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
        coordinates
            pd.DataFrame, N x 2 coordinate matrix. Rows correspond to spatial locations and columns correspond to the x and
            y spatial coordinates.

        Raises
        ----------
        Exception
            Raises excpetion if the number of spatial locations do not match between counts and coordinates.

    """
    
    if counts.shape[0] != coordinates.shape[0]:
        raise Exception("Number of spatial locations do not match between counts and coordinates!")
