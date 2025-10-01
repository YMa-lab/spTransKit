from typing import Literal

import numpy as np
import scanpy as sc
import pandas as pd

import rpy2.robjects as ro


def size_factor(matrix : np.ndarray) -> np.ndarray:
    """ This function calculates the library size factor for each spatial location.
    
    Parameters
    ----------
    matrix
        np.ndarray, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.

    Returns
    ----------
    size_factors
        np.ndarray, N x 1 array of library size factors (s) for each spatial location.

    """

    sf = matrix.sum(1) / (matrix.sum() / matrix.shape[0])

    return np.reshape(sf, (matrix.shape[0], 1))


def cpm(matrix : np.ndarray) -> np.ndarray:
    """ This function calculates the counts per million (CPM) library size factor for each spatial location.
    
    Parameters
    ----------
    matrix
        np.ndarray, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.

    Returns
    ----------
    cpm_size_factors
        np.ndarray, N x 1 array of CPM library size factors (s) for each spatial location.

    """

    return np.reshape((matrix.sum(1) / 1000000), (matrix.shape[0], 1))


def n(matrix : np.ndarray) -> np.ndarray:
    """ This function calculates the sum of all gene counts for each spatial location. Utilized for the Analytic Pearson
    transformations.
    
    Parameters
    ----------
    matrix
        np.ndarray, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.

    Returns
    ----------
    n
        np.ndarray, N x 1 array of gene count sums for each spatial location

    """

    return np.reshape(matrix.sum(1), (matrix.shape[0], 1))


def p(matrix : np.ndarray) -> np.ndarray:
    """ This function calculates the proportion of expression of each gene compared to the total gene expression across
    the entire tissue. Utilized for the Analytic Pearson transformations.
    
    Parameters
    ----------
    matrix
        np.ndarray, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.

    Returns
    ----------
    p
        np.ndarray, G x 1 array of expression proportions for each gene

    """

    return matrix.sum(0) / n(matrix).sum()


def handle_r_params(params : dict):
    """ This functions handles the conversion of parameters from Python data types to R data types using the rpy2
    package. Used in the SCTransform, Dino, and SpaNorm transformations.

    Parameters
    ----------
    params
        dict, dictionary of parameters passed into the parent function
    
    """

    for param in params.keys():
        if param == "matrix":
            continue
        if str(params[param]) == str(np.nan):
            params[param] = ro.r("NA")
        match params[param]:
            case None:
                params[param] = ro.r("NULL")
            case bool():
                if params[param] == True:
                    params[param] = ro.r("TRUE")
                else:
                    params[param] = ro.r("FALSE")
            case int():
                params[param] = ro.r(str(params[param]))
            case float():
                params[param] = ro.r(str(params[param]))
            case list():
                params[param] = ro.r.c(params[param])
            case dict():
                params[param] = ro.ListVector(params[param])
            case pd.DataFrame():
                continue
        if params[param] == "-Inf":
            params[param] = ro.r("-Inf")


def get_unfiltered_dlpfc_data(sample : Literal["151507", "151508", "151509", "151510", "151669", "151670", 
                                               "151671", "151672", "151673", "151674", "151675", "151676"]):
    """ This function retrieves the unfiltered human dorsolateral prefrontal cortex spatial transcriptomics data 
    gerenated by the 10x Visium platform (Maynard et al., 2021).
    
    Parameters
    ----------
    sample
        str, specifies the sample from the Human DLPFC dataset, options = 151507, 151508, 151509, 151510, 151669, 151670,
        151671, 151672, 151673, 151674, 151675, 151676

    Returns
    ----------
    matrix
        sc.AnnData, unflitered N x G gene count matrix. Rows correspond to spatial locations and columns correspond to 
        genes.
    coordinates
        pd.DataFrame, unflitered N x 2 coordinate matrix. Rows correspond to spatial locations and columns correspond to 
        the x and y spatial coordinates.
    
    """

    matrix = sc.read_10x_h5("data/DLPFC_" + sample + "/unfiltered_raw.h5")
    coordinates = pd.read_csv("data/DLPFC_" + sample + "/unfiltered_coord.txt", header = None, index_col = 0)

    coordinates = coordinates.loc[coordinates[1] != 0, 2:3]
    coordinates.columns = ["x", "y"]
    matrix = sc.AnnData(matrix.to_df().loc[coordinates.index, :])

    return matrix, coordinates


def get_filtered_dlpfc_data(sample : Literal["151507", "151508", "151509", "151510", "151669", "151670", 
                                             "151671", "151672", "151673", "151674", "151675", "151676"]):
    """ This function retrieves the filtered human dorsolateral prefrontal cortex spatial transcriptomics data generated
    by the 10x Visium platform (Maynard et al., 2021).
    
    Parameters
    ----------
    sample
        str, specifies the sample from the Human DLPFC dataset, options = 151507, 151508, 151509, 151510, 151669, 151670,
        151671, 151672, 151673, 151674, 151675, 151676

    Returns
    ----------
    matrix
        sc.AnnData, filtered N x G gene count matrix. Rows correspond to spatial locations and columns correspond to 
        genes.
    coordinates
        pd.DataFrame, fltered N x 2 coordinate matrix. Rows correspond to spatial locations and columns correspond to the 
        x and y spatial coordinates.
    
    """

    matrix = sc.read_h5ad("data/DLPFC_" + sample + "/filtered_raw.h5ad")
    coordinates = pd.read_csv("data/DLPFC_" + sample + "/filtered_coord.csv", header = 0, index_col = 0)

    return matrix, coordinates