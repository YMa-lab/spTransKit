import pandas as pd
import numpy as np
import scanpy as sc
import scanpy.preprocessing._deprecated as scppd
import normalisr.normalisr as n

import sptranskit.helpers as h

# Size Factor-Based Transformations

def size(matrix : sc.AnnData, inplace : bool = False):
    """ This function applies the y/s transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix: sc.AnnData, N x G gene count matrix
    inplace: bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy: sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """
    
    if inplace:
        matrix.X = np.divide(matrix.X, h.size_factor(matrix.X))
    else:
        copy = matrix.copy()
        copy.X = np.divide(matrix.X, h.size_factor(matrix.X))
        return copy

def cpm(matrix : sc.AnnData, inplace : bool = False):
    """ This function applies the CPM transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix: sc.AnnData, N x G gene count matrix
    inplace: bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy: sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    if inplace:
        matrix.X = np.divide(matrix.X, h.cpm(matrix.X))
    else:
        copy = matrix.copy()
        copy.X = np.divide(matrix.X, h.cpm(matrix.X))
        return copy
    
def weinreb(matrix : sc.AnnData, inplace : bool = False):
    """ This function applies the scanpy transformation to the gene count matrix according to Weinreb et al., 2017.
    
    Parameters
    ----------
    matrix: sc.AnnData, N x G gene count matrix
    inplace: bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy: sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    if inplace:
        sc.pp.log1p(matrix.X)
        matrix.X = scppd.normalize_per_cell_weinreb16_deprecated(matrix.X)
    else: 
        copy = matrix.copy()
        sc.pp.log1p(copy.X)
        copy.X = scppd.normalize_per_cell_weinreb16_deprecated(copy.X)
        return copy
    
def zheng(matrix : sc.AnnData, inplace : bool = False):
    """ This function applies the scanpy transformation to the gene count matrix according to Zheng et al., 2016.
    
    Parameters
    ----------
    matrix: sc.AnnData, N x G gene count matrix
    inplace: bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy: sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    if inplace:
        sc.pp.filter_genes(matrix, min_counts = 1)
        sc.pp._recipes.normalize_total(matrix, key_added = "n_counts_all")
        sc.pp.log1p(matrix.X)
        sc.pp.scale(matrix.X)
    else:
        copy = matrix.copy()
        sc.pp.filter_genes(copy, min_counts = 1)
        sc.pp._recipes.normalize_total(copy, key_added = "n_counts_all")
        sc.pp.log1p(copy.X)
        sc.pp.scale(copy.X)
        return copy
    
def tmm(matrix : sc.AnnData, trim_fraction : float = 0.2, inplace : bool = False):
    """ This function applies the TMM transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix: scanpy.AnnData, N x G gene count matrix
    trim_fraction: float, fraction of data to trim from the ends, default = 0.2
    inplace: bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy: sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    m = pd.DataFrame(matrix.X)

    ref_row = m.sum(axis = 1).sort_values().index.to_list()[m.shape[0] // 2]
    ref_counts = m.loc[ref_row, :]
    ref_counts = np.array(ref_counts / ref_counts.sum())
    ref_counts[ref_counts == 0.0] = np.nan
    norm_factors = [1] * len(m.index.to_list())
        
    for i in range(0, m.shape[0]):
        current_counts = m.iloc[i, :]
        current_counts = np.array(current_counts / current_counts.sum())
        current_counts[current_counts == 0] = np.nan
            
        log_ratios = np.log2((current_counts) / (ref_counts))
        abs_log_ratios = np.abs(log_ratios)
        weights = (1 / (current_counts)) + (1 / (ref_counts))
            
        finite = np.isfinite(log_ratios)
        log_ratios = log_ratios[finite]
        abs_log_ratios = abs_log_ratios[finite]
        weights = weights[finite]
            
        lower_trim = int(np.floor(trim_fraction * len(log_ratios)))
        upper_trim = int(np.ceil((1 - trim_fraction) * len(log_ratios)))
        sorted_indices = np.argsort(abs_log_ratios)
        trimmed_indices = sorted_indices[lower_trim:upper_trim]
        trimmed_log_ratios = log_ratios[trimmed_indices]
        trimmed_weights = weights[trimmed_indices]
        if np.sum(trimmed_weights) == 0.0:
            norm_factors[i] = 1.0
        else:
            norm_factors[i] = 2 ** (np.sum(trimmed_log_ratios * trimmed_weights) / np.sum(trimmed_weights))
        
    norm_factors = pd.Series(norm_factors, index = m.index)
    
    if inplace:
        matrix.X = np.array(m.div(norm_factors, axis = 0))
    else:
        copy = matrix.copy()
        copy.X = np.array(m.div(norm_factors, axis = 0))
        return copy
    
def deseq2(matrix : sc.AnnData, inplace : bool = False):
    """ This function applies the DESeq2 transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix: sc.AnnData, N x G gene count matrix
    inplace: bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy: sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    m = pd.DataFrame(matrix.X)

    log_counts = pd.DataFrame(np.log(m + 1), index = m.index, columns = m.columns)
    log_means = log_counts.mean(0)
    filtered_genes = log_means.loc[log_means > 0.5].index.to_list()
    log_ratios = log_counts.loc[:, filtered_genes].sub(log_means.loc[filtered_genes], axis = 1)

    log_medians = log_ratios.median(axis = 1)
    size_factors = pd.Series(np.exp(log_medians), index = log_medians.index)
    
    if inplace:
        matrix.X = np.array(m.div(size_factors, axis = 0))
    else:
        copy = matrix.copy()
        copy.X = np.array(m.div(size_factors, axis = 0))
        return copy
    

# Delta Method-Based Transformations

def shifted_log(matrix : sc.AnnData, inplace : bool = False):
    """ This function applies the log(y/s + 1) transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix: sc.AnnData, N x G gene count matrix
    inplace: bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy: sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    shift = np.divide(matrix.X, h.size_factor(matrix.X)) + 1
    
    if inplace:
        matrix.X = np.log10(shift)
    else:
        copy = matrix.copy()
        copy.X = np.log10(shift)
        return copy
    
def cpm_shifted_log(matrix : sc.AnnData, inplace : bool = False):
    """ This function applies the log(CPM + 1) transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix: sc.AnnData, N x G gene count matrix
    inplace: bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy: sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    shift = np.divide(matrix.X, h.cpm(matrix.X)) + 1
    
    if inplace:
        matrix.X = np.log10(shift)
    else:
        copy = matrix.copy()
        copy.X = np.log10(shift)
        return copy

def shifted_log_size(matrix : sc.AnnData, inplace : bool = False):
    """ This function applies the log(y/s + 1)/u transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix: sc.AnnData, N x G gene count matrix
    inplace: bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy: sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    u = h.size_normalization(matrix.X)
    shift = np.divide(matrix.X, h.size_factor(matrix.X)) + 1

    if inplace:
        matrix.X = np.divide(np.log10(shift), u)
    else:
        copy = matrix.copy()
        copy.X = np.divide(np.log10(shift), u)
        return copy

def acosh(matrix : sc.AnnData, alpha : float = 0.05, inplace : bool = False):
    """ This function applies the acosh(2αy/s + 1) transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix: sc.AnnData, N x G gene count matrix
    alpha: float, overdispersion hyperparamter, default = 0.05
    inplace: bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy: sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    size = np.divide(matrix.X, h.size_factor(matrix.X))
    shift = (size * alpha * 2) + 1

    if inplace:
        matrix.X = np.arccosh(shift)
    else:
        copy = matrix.copy()
        copy.X = np.arccosh(shift)
        return copy
    
def pseudo_shifted_log(matrix : sc.AnnData, alpha : float = 0.05, inplace : bool = False):
    """ This function applies the log(y/s + 1/4α) transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix: sc.AnnData, N x G gene count matrix
    alpha: float, overdispersion hyperparamter, default = 0.05
    inplace: bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy: sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    shift = np.divide(matrix.X, h.size_factor(matrix.X)) + (1 / (4 * alpha))

    if inplace:
        matrix.X = np.log10(shift)
    else:
        copy = matrix.copy()
        copy.X = np.log10(shift)
        return copy
    
# Model-Based Transformations

def analytic_pearson(matrix : sc.AnnData, alpha : float = 0.05, inplace : bool = False):
    """ This function applies the Analytic Pearson residual transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix: sc.AnnData, N x G gene count matrix
    alpha: float, overdispersion hyperparamter, default = 0.05
    inplace: bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy: sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    nb_matrix = np.zeros(matrix.X.shape)
    nb_matrix = np.add(nb_matrix, h.n(matrix.X))
    nb_matrix = np.multiply(nb_matrix, h.p(matrix.X))

    var_matrix = nb_matrix.copy()
    var_matrix = np.add(((var_matrix ** 2) * alpha), nb_matrix) ** 0.5

    final_matrix = np.divide((np.subtract(matrix.X, nb_matrix)), var_matrix)

    if inplace:
        matrix.X = final_matrix
    else:
        copy = matrix.copy()
        copy.X = final_matrix
        return copy

def sc_pearson(matrix : sc.AnnData, inplace : bool = False):
    """ This function applies the scanpy Pearson residual transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix: sc.AnnData, N x G gene count matrix
    inplace: bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy: sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    if inplace:
        sc.experimental.pp.normalize_pearson_residuals(matrix)
    else:
        copy = matrix.copy()
        sc.experimental.pp.normalize_pearson_residuals(copy)
        return copy
    
def seurat(matrix : sc.AnnData, inplace : bool = False):
    """ This function applies the Seurat transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix: sc.AnnData, N x G gene count matrix
    inplace: bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy: sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    if inplace:
        sc.pp.normalize_total(matrix, target_sum = 1e4)
        sc.pp.log1p(matrix)
        sc.pp.scale(matrix, max_value = 10)
    else:
        copy = matrix.copy()
        sc.pp.normalize_total(copy, target_sum = 1e4)
        sc.pp.log1p(copy)
        sc.pp.scale(copy, max_value = 10)
        return copy

def normalisr(matrix : sc.AnnData, inplace : bool = False):
    """ This function applies the Normalisr transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix: sc.AnnData, N x G gene count matrix
    inplace: bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy: sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    norm = n.lcpm(np.transpose(matrix.X))[0]

    if inplace:
        matrix.X = np.transpose(norm)
    else:
        copy = matrix.copy()
        copy.X = np.transpose(norm)
        return copy
    
def psinorm(matrix : sc.AnnData, inplace : bool = False):
    """ This function applies the PsiNorm transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix: sc.AnnData, N x G gene count matrix
    inplace: bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy: sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    n = matrix.shape[0]
    m = (matrix.X + 1).min(0)
    log_counts = np.log(matrix.X + 1)
    log_m = np.log(m)
    sf = np.subtract(log_counts, log_m).sum(1) / n
    sf = np.reshape(sf, (matrix.shape[0], 1))
    
    if inplace:
        matrix.X = np.divide(matrix.X, sf)
    else:
        copy = matrix.copy()
        copy.X = np.divide(matrix.X, sf)
        return copy
