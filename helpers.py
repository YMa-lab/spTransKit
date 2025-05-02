import numpy as np

def size_factor(matrix : np.ndarray) -> np.ndarray:
    """ This function calculates the library size factor for each spatial location.
    
    Parameters
    ----------
    matrix: np.ndarray, N x G gene count matrix

    Returns
    ----------
    size_factors: np.ndarray, N x 1 array of library size factors (s) for each spatial location

    """

    return np.reshape(np.divide(matrix.sum(1), (matrix.sum() / matrix.shape[0])), (matrix.shape[0], 1))

def cpm(matrix : np.ndarray) -> np.ndarray:
    """ This function calculates the counts per million (CPM) library size factor for each spatial location.
    
    Parameters
    ----------
    matrix: np.ndarray, N x G gene count matrix

    Returns
    ----------
    cpm_size_factors: np.ndarray, N x 1 array of CPM library size factors (s) for each spatial location

    """

    return np.reshape((matrix.sum(1) / 1000000), (matrix.shape[0], 1))

def size_normalization(matrix : np.ndarray) -> np.ndarray:
    """ This function calculates the sum of all gene counts for each spatial location. Utilized for the log(y/s + 1)/u
    transformation.
    
    Parameters
    ----------
    matrix: np.ndarray, N x G gene count matrix

    Returns
    ----------
    u: np.ndarray, N x 1 array of gene count sums for each spatial location

    """

    return np.reshape(matrix.sum(1), (matrix.shape[0], 1))

def n(matrix : np.ndarray) -> np.ndarray:
    """ This function calculates the sum of all gene counts for each spatial location. Utilized for the Analytic Pearson
    transformation.
    
    Parameters
    ----------
    matrix: np.ndarray, N x G gene count matrix

    Returns
    ----------
    n: np.ndarray, N x 1 array of gene count sums for each spatial location

    """

    return np.reshape(matrix.sum(1), (matrix.shape[0], 1))

def p(matrix : np.ndarray) -> np.ndarray:
    """ This function calculates the proportion of expression of each gene compared to the total gene expression across
    the entire tissue. Utilized for the Analytic Pearson transformation.
    
    Parameters
    ----------
    matrix: np.ndarray, N x G gene count matrix

    Returns
    ----------
    p: np.ndarray, G x 1 array of expression proportions for each gene

    """

    return np.divide(matrix.sum(0), n(matrix).sum())
