from typing import Literal

import pandas as pd
import numpy as np
import scanpy as sc
import normalisr.normalisr as n

import rpy2.robjects as ro
import rpy2.robjects.packages as rp
import rpy2.robjects.vectors as rv
import rpy2.robjects.pandas2ri as rpd
import rpy2.robjects.numpy2ri as rnp
base = rp.importr("base")
utils = rp.importr("utils")

import sptranskit.helpers as h


# Size Factor-Based Transformations

def size(data : sc.AnnData, inplace : bool = False):
    """ This function applies the y/s transformation to the gene count matrix.
    
    Parameters
    ----------
    data
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
        Spatial information should be stored in data.obsm["spatial"]
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False.

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False.

    """

    if data.obsm["spatial"] is None:
        raise Exception("No spatial information included in AnnData object.")

    data.raw = data.copy()
    
    if inplace:
        data.X = np.divide(data.X, h.size_factor(data.X))
    else:
        copy = data.copy()
        copy.X = np.divide(data.X, h.size_factor(data.X))
        return copy


def cpm(data : sc.AnnData, inplace : bool = False):
    """ This function applies the CPM transformation to the gene count matrix.
    
    Parameters
    ----------
    data
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
        Spatial information should be stored in data.obsm["spatial"]
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False.

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False.

    """

    if data.obsm["spatial"] is None:
        raise Exception("No spatial information included in AnnData object.")

    data.raw = data.copy()

    if inplace:
        data.X = np.divide(data.X, h.cpm(data.X))
    else:
        copy = data.copy()
        copy.X = np.divide(data.X, h.cpm(data.X))
        return copy
    

def zheng(data : sc.AnnData, inplace : bool = False):
    """ This function applies the scanpy transformation to the gene count matrix according to Zheng et al., 2016.
    
    Parameters
    ----------
    data
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
        Spatial information should be stored in data.obsm["spatial"]
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False.

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False.

    """

    if data.obsm["spatial"] is None:
        raise Exception("No spatial information included in AnnData object.")

    data.raw = data.copy()

    if inplace:
        sc.pp.normalize_total(data, key_added = "n_counts_all")
        sc.pp.log1p(data.X)
        sc.pp.scale(data.X)
    else:
        copy = data.copy()
        sc.pp.normalize_total(copy, key_added = "n_counts_all")
        sc.pp.log1p(copy.X)
        sc.pp.scale(copy.X)
        return copy
    

def tmm(data : sc.AnnData,
        logratio_trim : float = 0.3,
        sum_trim : float = 0.05,
        inplace : bool = False):
    """ This function applies the TMM transformation to the gene count matrix.
    
    Parameters
    ----------
    data
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
        Spatial information should be stored in data.obsm["spatial"]
    logratio_trim
        float, fraction of data to trim from the ends of ranked log ratios, default = 0.3.
    sum_trim
        float, fraction of data to trim from the ends of ranked sums, default = 0.05.
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False.

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False.

    """

    if data.obsm["spatial"] is None:
        raise Exception("No spatial information included in AnnData object.")

    np.seterr(all = "ignore")

    m = pd.DataFrame(data.X)

    f = m.quantile(q = 0.75, axis = 1).div(m.sum(1))
    ref_index = (abs(f - f.mean())).argmin()
    ref_row = m.iloc[ref_index, :]

    norm_factors = pd.Series(1.0, index = m.index)
        
    for i in range(0, len(norm_factors.index)):
        current_row = m.iloc[i, :]
            
        log_ratios = np.log2((current_row / current_row.sum()) / (ref_row / ref_row.sum()))
        abs_exp = (np.log2(current_row / current_row.sum()) + np.log2(ref_row / ref_row.sum())) / 2
        weights = ((current_row.sum() - current_row) / current_row.sum() / current_row) + ((ref_row.sum() - ref_row) / ref_row.sum() / ref_row)
            
        keep = np.logical_and(np.logical_and(np.isfinite(log_ratios), np.isfinite(abs_exp)), (abs_exp > -1e10))
        log_ratios = log_ratios[keep]
        abs_exp = abs_exp[keep]
        weights = weights[keep]

        if log_ratios.max() < 1e-6:
            continue

        n = len(log_ratios)
        low_l = np.floor(n * logratio_trim) + 1
        high_l = n + 1 - low_l
        low_s = np.floor(n * sum_trim) + 1
        high_s = n + 1 - low_s
            
        keep_l = np.logical_and(log_ratios.rank() >= low_l, log_ratios.rank() <= high_l)
        keep_s = np.logical_and(abs_exp.rank() >= low_s, abs_exp.rank() <= high_s)
        keep = np.logical_and(keep_l, keep_s)
        
        log_ratios = log_ratios[keep]
        weights = weights[keep]

        sf = (log_ratios / weights).sum(skipna = True) / (1 / weights).sum(skipna = True)
        
        if str(sf) == str(np.nan):
            continue
        else:
            norm_factors.iloc[i] = 2 ** sf
        
    norm_factors = pd.Series(norm_factors, index = m.index)

    data.raw = data.copy()
    
    if inplace:
        data.X = np.array(m.div(norm_factors, axis = 0))
    else:
        copy = data.copy()
        copy.X = np.array(m.div(norm_factors, axis = 0))
        return copy
    

def deseq2(data : sc.AnnData, inplace : bool = False):
    """ This function applies the DESeq2 transformation to the gene count matrix.
    
    Parameters
    ----------
    data
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
        Spatial information should be stored in data.obsm["spatial"]
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False.

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False.

    """

    if data.obsm["spatial"] is None:
        raise Exception("No spatial information included in AnnData object.")

    np.seterr(all = "ignore")

    m = pd.DataFrame(data.X)

    log_counts = pd.DataFrame(np.log(m), index = m.index, columns = m.columns)
    log_counts.replace(float("-inf"), np.nan, inplace = True)
    log_means = log_counts.mean(0, skipna = True)
    filtered_genes = log_means.loc[log_means > 0].index.to_list()
    log_ratios = log_counts.loc[:, filtered_genes].sub(log_means.loc[filtered_genes], axis = 1)

    log_medians = log_ratios.median(axis = 1)
    log_medians.fillna(0, inplace = True)
    size_factors = pd.Series(np.exp(log_medians), index = log_medians.index)
    
    data.raw = data.copy()

    if inplace:
        data.X = np.array(m.div(size_factors, axis = 0))
    else:
        copy = data.copy()
        copy.X = np.array(m.div(size_factors, axis = 0))
        return copy
    

# Delta Method-Based Transformations

def shifted_log(data : sc.AnnData, inplace : bool = False):
    """ This function applies the log(y/s + 1) transformation to the gene count matrix.
    
    Parameters
    ----------
    data
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
        Spatial information should be stored in data.obsm["spatial"]
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False.

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False.

    """

    if data.obsm["spatial"] is None:
        raise Exception("No spatial information included in AnnData object.")
    
    shift = np.log1p((data.X / h.size_factor(data.X)))

    data.raw = data.copy()

    if inplace:
        data.X = shift
    else:
        copy = data.copy()
        copy.X = shift
        return copy
    

def cpm_shifted_log(data : sc.AnnData, inplace : bool = False):
    """ This function applies the log(CPM + 1) transformation to the gene count matrix.
    
    Parameters
    ----------
    data
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
        Spatial information should be stored in data.obsm["spatial"]
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False.

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False.

    """

    if data.obsm["spatial"] is None:
        raise Exception("No spatial information included in AnnData object.")

    shift = np.log1p(data.X / h.cpm(data.X))
    
    data.raw = data.copy()

    if inplace:
        data.X = shift
    else:
        copy = data.copy()
        copy.X = shift
        return copy


def shifted_log_size(data : sc.AnnData, inplace : bool = False):
    """ This function applies the log(y/s + 1)/u transformation to the gene count matrix.
    
    Parameters
    ----------
    data
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
        Spatial information should be stored in data.obsm["spatial"]
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False.

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False.

    """

    if data.obsm["spatial"] is None:
        raise Exception("No spatial information included in AnnData object.")

    shift = np.log1p(data.X / h.size_factor(data.X))

    data.raw = data.copy()

    if inplace:
        data.X = shift / h.size_factor(shift)
    else:
        copy = data.copy()
        copy.X = shift / h.size_factor(shift)
        return copy


def acosh(data : sc.AnnData, alpha : float = 0.05, inplace : bool = False):
    """ This function applies the acosh(2αy/s + 1) transformation to the gene count matrix.
    
    Parameters
    ----------
    data
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
        Spatial information should be stored in data.obsm["spatial"]
    alpha
        float, overdispersion hyperparamter, default = 0.05.
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False.

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False.

    """

    if data.obsm["spatial"] is None:
        raise Exception("No spatial information included in AnnData object.")

    size = data.X / h.size_factor(data.X)
    shift = (1 / np.sqrt(alpha)) * np.arccosh((size * alpha * 2) + 1)

    data.raw = data.copy()

    if inplace:
        data.X = shift
    else:
        copy = data.copy()
        copy.X = shift
        return copy
    

def pseudo_shifted_log(data : sc.AnnData, alpha : float = 0.05, inplace : bool = False):
    """ This function applies the log(y/s + 1/4α) transformation to the gene count matrix.
    
    Parameters
    ----------
    data
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
        Spatial information should be stored in data.obsm["spatial"]
    alpha
        float, overdispersion hyperparamter, default = 0.05.
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False.

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False.

    """

    if data.obsm["spatial"] is None:
        raise Exception("No spatial information included in AnnData object.")

    shift = (1 / np.sqrt(alpha)) * np.log1p(4 * alpha * (data.X / h.size_factor(data.X)))

    data.raw = data.copy()

    if inplace:
        data.X = shift
    else:
        copy = data.copy()
        copy.X = shift
        return copy
    

# Model-Based Transformations

def analytic_pearson_noclip(data : sc.AnnData, alpha : float = 0.05, inplace : bool = False):
    """ This function applies the Analytic Pearson (no clip) transformation to the gene count matrix.
    
    Parameters
    ----------
    data
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
        Spatial information should be stored in data.obsm["spatial"]
    alpha
        float, overdispersion hyperparamter, default = 0.05.
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False.

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False.

    """

    if data.obsm["spatial"] is None:
        raise Exception("No spatial information included in AnnData object.")

    nb_matrix = np.zeros(data.X.shape)
    nb_matrix = np.add(nb_matrix, h.n(data.X))
    nb_matrix = np.multiply(nb_matrix, h.p(data.X))

    var_matrix = nb_matrix.copy()
    var_matrix = np.add(((var_matrix ** 2) * alpha), nb_matrix) ** 0.5

    final_matrix = np.divide((np.subtract(data.X, nb_matrix)), var_matrix)

    data.raw = data.copy()

    if inplace:
        data.X = final_matrix
    else:
        copy = data.copy()
        copy.X = final_matrix
        return copy


def analytic_pearson_clip(data : sc.AnnData, alpha : float = 0.05, clip : float = None, inplace : bool = False):
    """ This function applies the Analytic Pearson (clip) transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
        Spatial information should be stored in data.obsm["spatial"]
    alpha
        float, overdispersion hyperparamter, default = 0.05.
    clip
        float, handles how the Pearson residuals are clipped. 
        • If None, resiudals are clipped using the default interval [-√N, √N], where N is the number of spatial locations
          in the dataset.
        • If a scalar, $c$ is provided, then the residuals will be clipped using the interval [-c, c].
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False.

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False.

    """

    if data.obsm["spatial"] is None:
        raise Exception("No spatial information included in AnnData object.")

    nb_matrix = np.zeros(data.X.shape)
    nb_matrix = np.add(nb_matrix, h.n(data.X))
    nb_matrix = np.multiply(nb_matrix, h.p(data.X))

    var_matrix = nb_matrix.copy()
    var_matrix = np.add(((var_matrix ** 2) * alpha), nb_matrix) ** 0.5

    final_matrix = np.divide((np.subtract(data.X, nb_matrix)), var_matrix)
    
    if clip is None:
        clip = matrix.X.shape[0] ** 0.5
    
    final_matrix = np.where(final_matrix > clip, clip, final_matrix)
    final_matrix = np.where(final_matrix < -clip, -clip, final_matrix)

    data.raw = data.copy()
    
    if inplace:
        data.X = final_matrix
    else:
        copy = data.copy()
        copy.X = final_matrix
        return copy


def sc_pearson(data : sc.AnnData, 
               theta : float = 100,
               clip : float = None,
               check_values : bool = True,
               layer : str = None,
               inplace : bool = False):
    """ This function applies the scanpy Pearson residual transformation to the gene count matrix.
    
    Parameters
    ----------
    data
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
        Spatial information should be stored in data.obsm["spatial"]
    theta
        float, the negative binomial overdispersion parameter for Pearson residuals, equal to 1/α, the overdisperson
        parameter used in the other Pearson residual transformations.
    clip
        float, determines if and how residuals should be clipped. 
        • If None, residuals are clipped to the interval [-√N, √N], where N is the number of spatial locations in the
          dataset.
        • If any scalar, c, is provided, residuals are clipped to the interval [-c, c]. Set clip = np.inf for no
          clipping.
    check_values
        bool, checks to see if counts are integers, default = True.
        • If True, checks if counts in selected layer are integers as expected by this function, and return a warning if
          non-integers are found.
        • If False, proceeds without checking. This can speed up code for large datasets.
    layer
        str, layer to use as input instead of X, default = None, in which X is used.
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False.

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False.

    """

    if data.obsm["spatial"] is None:
        raise Exception("No spatial information included in AnnData object.")

    data.raw = data.copy()

    if inplace:
        sc.experimental.pp.normalize_pearson_residuals(data, theta = theta, clip = clip, check_values = check_values, 
                                                       layer = layer)
    else:
        copy = data.copy()
        sc.experimental.pp.normalize_pearson_residuals(copy)
        return copy


def normalisr(data : sc.AnnData, 
              normalize : bool = True,
              nth : int = 0,
              ntot : int = None,
              varscale : float = 0,
              seed : int = None,
              lowmem : bool = True,
              nocov : bool = False,
              inplace : bool = False):
    """ This function applies the Normalisr transformation to the gene count matrix.
    
    Parameters
    ----------
    data
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
        Spatial information should be stored in data.obsm["spatial"]
    normalize,
        bool, whether to normalize output to logCPM per cell, defualt = True.
    nth
        int, number of threads to use, default = 0, which uses all cores automatically detected.
    ntot
        int, manually sets value of total number of reads in binomial distribution, default = None. Since the posterior 
        distribution stabilizes quickly as ntot increases, a large number (e.g. 1e9) is good for general use. The default
        of None disables the manual value.
    varscale
        float, resamples estimated expression using the posterior Beta distribution, default = 0. varscale sets the scale
        of variance rather than its actual value from the posterior distribution.
    seed
        int, inital random seed if set, default = None.
    lowmem
        bool, low memory mode, default = True.
    nocov
        bool, whether to skip producing covariate variables, default = False.
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False.

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False.

    """

    if data.obsm["spatial"] is None:
        raise Exception("No spatial information included in AnnData object.")

    norm = n.lcpm(np.transpose(data.X), normalize = normalize, nth = nth, ntot = ntot, varscale = varscale, seed = seed,
                  lowmem = lowmem, nocov = nocov)[0]
    
    data.raw = data.copy()

    if inplace:
        data.X = np.transpose(norm)
    else:
        copy = data.copy()
        copy.X = np.transpose(norm)
        return copy
    

def psinorm(data : sc.AnnData, inplace : bool = False):
    """ This function applies the PsiNorm transformation to the gene count matrix.
    
    Parameters
    ----------
    data
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
        Spatial information should be stored in data.obsm["spatial"]
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False.

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False.

    """

    if data.obsm["spatial"] is None:
        raise Exception("No spatial information included in AnnData object.")

    n = data.shape[0]
    m = (data.X + 1).min(1)
    log_counts = np.log(data.X + 1)
    log_m = np.log(m)
    log_m = np.reshape(log_m, (data.shape[0], 1))
    sf = np.subtract(log_counts, log_m).sum(1) / n
    sf = np.reshape(sf, (data.shape[0], 1))

    data.raw = data.copy()
    
    if inplace:
        data.X = np.divide(data.X, sf)
    else:
        copy = data.copy()
        copy.X = np.divide(copy.X, sf)
        return copy


def sctransform(data : sc.AnnData, 
                cell_attr : pd.DataFrame = None,
                latent_var : list[str] = ["log_umi"],
                batch_var : list[str] = None,
                latent_var_nonreg : list[str] = None,
                n_genes : int = 2000,
                n_cells : int = None,
                method : Literal["poisson", "qpoisson", "nb_fast", "nb", "nb_theta_given", "glmGamPoi", "offset",
                                 "offset_shared_theta_estimate", "glmGamPoi_offset"] = "poisson",
                do_regularize : bool = True,
                theta_regularization : Literal["od_factor", "log_theta"] = "od_factor",
                res_clip_range : list[float] = None,
                bin_size : int = 500,
                min_cells : int = 5,
                residual_type : Literal["pearson", "deviance", "none"] = "pearson",
                return_cell_attr : bool = False,
                return_gene_attr : bool = True,
                return_corrected_umi : bool = False,
                min_variance : float | Literal["umi_median", "model_median", "model_mean", "-Inf"] = "-Inf",
                bw_adjust : float = 3,
                gmean_eps : float = 1,
                theta_estimation_fun : Literal["theta.ml", "theta.mm"] = "theta.ml",
                theta_given : list | float = None,
                exclude_poisson : bool = False,
                use_geometric_mean : bool = True,
                use_geometric_mean_offset : bool = False,
                fix_intercept : bool = False,
                fix_slope : bool = False,
                scale_factor : float = np.nan,
                vst_flavor : Literal["v2"] = None,
                verbosity : int = 2,
                inplace : bool = False):
    """ This function applies the SCTransform transformation to the gene count matrix.
    
    Parameters
    ----------
    data
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
        Spatial information should be stored in data.obsm["spatial"]
    cell_attr
        pd.DataFrame, a data frame containining the dependent variables, default = None, in which a data frame with UMI
        and genes is generated.
    latent_var
        list[str], the independent variables to regress out, expressed as a list of strings, default = ["log_umi"]. These
        names must match the column names in cell_attr.
    batch_var
        list[str], the independent variables indicating which batch a spatial location belongs to, default = None, in
        which no batch interaction terms are used.
    latent_var_nonreg
        list[str], the non-regularized dependent variables to regress out, expressed as a list of strings, default = 
        None. These names must match the column names in cell_attr.
    n_genes
        int, number of genes to use when estimating parameters, default = 2000.
    n_cells
        int, number of spatial locations to use when estimating parameters, default = None, in which all spatial
        locations are used.
    method
        str, method to use for inital parameter estimation. Can be one of "poisson", "qpoisson", "nb_fast", "nb",
        "nb_theta_given", "glmGamPoi", "offset", "offset_shared_theta_estimate", or "glmGamPoi_offset", default = 
        "poisson".
    do_regularize
        bool, whether to use the n_genes parameter in initial parameter estimation, default = True.
    theta_regularization
        str, method to use to regularize theta, default = "od_factor"
    res_clip_range
        list[float], list of numbers of length 2, specifying the min and max values the results will be clipped to,
        default = None, which uses the default interval [-√N, √N], where N is the number of spatial locations in the
        dataset.
    bin_size
        int, number of genes to process simultaneously, default = 500.
    min_cells
        int, only use genes that have been detected in this many cells, default = 5.
    residual_type
        str, which type of residuals to return. Can be one of "pearson", "deviance", or "none", default = "pearson".
    return_cell_attr
        bool, make cell attributes part of the output, default = False.
    return_gene_attr
        bool, calculate gene attributes and make part of the output, default = True.
    return_corrected_umi
        bool, return output containing corrected UMI matric, defualt = False.
    min_variance
        float or str, lower bound for the estimated variance for any gene in any cell when calculating the Pearson
        residual. Can be one of "umi_median", "model_median", "model_mean", "-Inf", or a scalar, default = "-Inf".
        • If "umi_median", uses (median of non-zero UMIs / 5)^2 as the minimum variance.
        • If "model_median" or "model_mean", uses the median/mean of the model estimated µ per gene as the minimum
          variance.
    bw_adjust
        float, kernel bandwidth adjustment factor used during regularization, default = 3.
    gmean_eps
        float, small value added when calculating geometric mean of a gene to avoid log(0), default = 1.
    theta_estimation_fun
        str, indicates which method to use to estimate theta when method = poisson, default = "theta.ml".
    theta_given
        list or float, specified values for theta parameter, default = None.
        • If method = "nb_theta_given", should be a list of numbers of fixed theta values for each gene.
        • If method = "offset", should be a single value.
    exclude_poisson
        bool, exclude poisson genes (e.g. µ < 0.001 or µ > variance) from regularization, default = False.
    use_geometric_mean
        bool, use geometric mean instead of arithmetic mean for all calculations, default = True.
    use_geometric_mean_offset
        bool, use geometric mean instead of arithmetic mean in the offset model for all calculations, default = False.
    fix_intercept
        bool, fix intercept as defined in the offset model, default = False.
    fix_slope
        bool, fix slope to log(10), default = False.
    scale_factor
        float, replace all values of UMI in the regression model by this value instead of the median UMI, default = 
        np.nan.
    vst_flavor
        str, flavor of the model to use, default = None.
        • If None, which uses the original sctransform model
        • If "v2", sets method = "glmGamPoi_offset", n_cells = 2000, and exclude_poisson = True.
    verbosity
        int, sets message output settings, default = 2.
        • If 0, shows no messages when function is running.
        • If 1, shows only messages when function is running.
        • If 2, shows both messages and progress bars when function is running.
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False.

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False.

    """

    if data.obsm["spatial"] is None:
        raise Exception("No spatial information included in AnnData object.")

    packages = ["sctransform"]
    names_to_install = [x for x in packages if not rp.isinstalled(x)]
    if len(names_to_install) > 0:
        utils.install_packages(rv.StrVector(names_to_install))
    for x in packages:
        rp.importr(x)

    params = locals()
    h.handle_r_params(params)
    if res_clip_range is None:
        params["res_clip_range"] = ro.FloatVector([-(data.shape[0] ** 0.5), data.shape[0] ** 0.5])
    else:
        params["res_clip_range"] = ro.FloatVector(res_clip_range)
    r_code = f"""function(x, cell_attr, latent_var, batch_var, latent_var_nonreg, n_genes, n_cells, method, 
    do_regularize, theta_regularization, res_clip_range, bin_size, min_cells, residual_type, return_cell_attr,
    return_gene_attr, return_corrected_umi, min_variance, bw_adjust, gmean_eps, theta_estimation_fun, theta_given,
    exclude_poisson, use_geometric_mean, use_geometric_mean_offset, fix_intercept, fix_slope, scale_factor, vst.flavor,
    verbosity) {{
        as.data.frame(vst(as.matrix(x),
            cell_attr = cell_attr,
            latent_var = latent_var,
            batch_var = batch_var,
            latent_var_nonreg = latent_var_nonreg,
            n_genes = n_genes,
            n_cells = n_cells,
            method = method,
            do_regularize = do_regularize,
            theta_regularization = theta_regularization,
            res_clip_range = res_clip_range,
            bin_size = bin_size,
            min_cells = min_cells,
            residual_type = residual_type,
            return_cell_attr = return_cell_attr,
            return_gene_attr = return_gene_attr,
            return_corrected_umi = return_corrected_umi,
            min_variance = min_variance,
            bw_adjust = bw_adjust,
            gmean_eps = gmean_eps,
            theta_estimation_fun = theta_estimation_fun,
            theta_given = theta_given,
            exclude_poisson = exclude_poisson,
            use_geometric_mean = use_geometric_mean,
            use_geometric_mean_offset = use_geometric_mean_offset,
            fix_intercept = fix_intercept,
            fix_slope = fix_slope,
            scale_factor = scale_factor,
            vst.flavor = vst.flavor,
            verbosity = verbosity)[["y"]])
    }}"""
    vst = ro.r(r_code)

    counts = data.to_df().T
    with (ro.default_converter + rpd.converter).context():
        trans = vst(counts, params["cell_attr"], params["latent_var"], params["batch_var"], params["latent_var_nonreg"],
                    params["n_genes"], params["n_cells"], params["method"], params["do_regularize"],
                    params["theta_regularization"], params["res_clip_range"], params["bin_size"], params["min_cells"],
                    params["residual_type"], params["return_cell_attr"], params["return_gene_attr"],
                    params["return_corrected_umi"], params["min_variance"], params["bw_adjust"], params["gmean_eps"],
                    params["theta_estimation_fun"], params["theta_given"], params["exclude_poisson"],
                    params["use_geometric_mean"], params["use_geometric_mean_offset"], params["fix_intercept"],
                    params["fix_slope"], params["scale_factor"], params["vst_flavor"], params["verbosity"])
    
    data.raw = data.copy()

    if inplace:
        data = data[trans.T.index.to_list(), trans.T.columns.to_list()]
        data.X = np.array(trans.T)
    else:
        copy = data.copy()
        copy = copy[trans.T.index.to_list(), trans.T.columns.to_list()]
        copy.X = np.array(trans.T)
        return copy


def dino(data : sc.AnnData, 
         nCores : int = 2,
         prec : int = 3,
         minNZ : int = 10,
         nSubGene : int = 10000,
         nSubCell : int = 10000,
         depth : list[float] = None,
         slope : float = None,
         minSlope : float = 0.5,
         maxSlope : float = 2,
         clusterSlope : bool = True,
         returnMeta : bool = False,
         doRQS : bool = False,
         emPar : dict = {"maxIter" : 100, "tol" : 0.1, "conPar" : 15, "maxK" : 100},
         inplace : bool = False):
    """ This function applies the Dino transformation to the gene count matrix.
    
    Parameters
    ----------
    data
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
        Spatial information should be stored in data.obsm["spatial"]
    nCores
        int, non-negative integer scalar denoting the number of cores to be used, default = 2. Setting nCores to 0, uses
        all cores that are automatically detected.
    prec
        int, positive integer denoting the number of decimals to which to round depth and normalized counts for
        computational efficiency, default = 3.
    minNZ
        int, positive integer denoting the minimum number of non-zero counts for a gene to be normalized by the Dino
        algorithm, default = 10.
    nSubGene
        int, positive number denoting the number of genes to subset for calculation of slope, default = 10000.
    nSubCell
        int, positive number denoting the number of spatial locations to subset for calculation of slope and the EM
        algorithm, default = 10000.
    depth
        list[float], list of length equal to number of spatial locations of median-centered, log-scaled measures of 
        cell-wise sequencing depth, default = None.
    slope
        float, scalar denoting the count-depth relationship on the log-log scale, default = None, which allows the Dino
        algorithm to calculate slope internally.
    minSlope
        float, scalar denoting the minimum slope, default = 0.5.
    maxSlope
        float, scalar denoting the maximum slope, defualt = 2.
    clusterSlope
        bool, whether cells should be pre-clustered prior to calculation of slope and cluster should then be used as a
        factor in regression, default = True.
    returnMeta
        bool, whether metadata (sequencing depth and slope) should be returned, default = False.
    doRQS
        bool, which method normalization resampling should be done, default = False.
        • If False, normalization is done by resampling the entire posterior distribution.
        • If True, restricted quantile sampling (RQS) can be performed to enforce stronger preservation of expression
          ranks in normalized data.
    emPar
        dict, parameters to send to EM algorithm
        • maxIter denotes the maximum number of model updates.
        • tol denotes the cutoff threshold for reductions in the log likelihood function.
        • conPar denotes the concentration parameter for the resampling.
        • maxK denotes the maximum number of mixture components in the mixture model.
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False.

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False.

    """

    if data.obsm["spatial"] is None:
        raise Exception("No spatial information included in AnnData object.")

    packages = ["Dino"]
    names_to_install = [x for x in packages if not rp.isinstalled(x)]
    if len(names_to_install) > 0:
        utils.install_packages(rv.StrVector(names_to_install))
    for x in packages:
        rp.importr(x)

    params = locals()
    h.handle_r_params(params)
    r_code = f"""function(x, nCores, prec, minNZ, nSubGene, nSubCell, depth, slope, minSlope, maxSlope, clusterSlope, 
    returnMeta, doRQS, emPar) {{
        as.data.frame(as.matrix(Dino(as.matrix(x),
                                     nCores = nCores,
                                     prec = prec, minNZ = minNZ,
                                     nSubGene = nSubGene,
                                     nSubCell = nSubCell,
                                     depth = depth,
                                     slope = slope,
                                     minSlope = minSlope,
                                     maxSlope = maxSlope,
                                     clusterSlope = clusterSlope,
                                     returnMeta = returnMeta,
                                     doRQS = doRQS,
                                     emPar = emPar)))
    }}"""
    dino = ro.r(r_code)

    counts = data.to_df().T
    with (ro.default_converter + rpd.converter).context():
        trans = dino(counts, params["nCores"], params["prec"], params["minNZ"], params["nSubGene"], params["nSubCell"], 
                     params["depth"], params["slope"], params["minSlope"], params["maxSlope"], params["clusterSlope"],
                     params["returnMeta"], params["doRQS"], params["emPar"])

    data.raw = data.copy()

    if inplace:
        data.X = trans.T
    else:
        copy = data.copy()
        copy.X = trans.T
        return copy


# Spatially Aware Transformations

def spanorm(data : sc.AnnData,
            sample_p : float = 0.25,
            gene_model : Literal["nb"] = "nb",
            adj_method : Literal["auto", "logpac", "pearson", "medbio", "meanbio"] = "auto",
            scale_factor : float = 1,
            df_tps : float = 6,
            lambda_a : float = 0.0001,
            batch : pd.DataFrame = None,
            tol : float = 0.0001,
            step_factor : float = 0.5,
            maxit_nb : int = 50,
            maxit_psi : int = 25,
            maxn_psi : int = 500,
            overwrite : bool = False,
            verbose : bool = True,
            inplace : bool = False):
    """ This function applies the SpaNorm transformation to the gene count matrix.
    
    Parameters
    ----------
    data
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
        Spatial information should be stored in data.obsm["spatial"]
    coordinates
        pd.DataFrame, N x 2 coordinate matrix. Rows correspond to spatial locations and columns correspond to the x and
        y spatial coordinates.
    sample_p
        float, maximum proportion of spatial locations to sample for model fitting, default = 0.25.
    gene_model
        str, model to use for gene/protein abundances, default = "nb".
    adj_method
        str, method to use to adjust the data. Can be one of "auto", "logpac", "pearson", "medbio", or "meanbio", default
        = "auto".
    scale_factor
        float, sample-specific scaling factor to adjust the counts, default = 1.
    df_tps
        float, degrees of freedom for the thin-plate spline, default = 6.
    lambda_a
        float, smoothing parameter for regularizing regression coefficients, default = 0.0001
    batch
        pd.DataFrame, specifies batch design to regress out, default = None.
    tol
        float, tolerance for convergence, default = 0.0001
    step_factor
        float, multiplicative factor to decrease IRLS step when log-likelihood diverges, default = 0.5.
    maxit_nb
        int, maximum number of IRLS iterations for estimating NB mean parameters for a given dispersion parameter,
        default = 50.
    maxit_psi
        int, maximum number of IRLS iterations for estimating the dispersion parameter, default = 25.
    maxn_psi
        int, maximum number of spatial locations to sample for dispersion estimation, default = 500.
    overwrite
        bool, whether to force recomputation and overwrite an existing fit, default = False.
    verbose
        bool, whether to show update messages, default = True.
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False.

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False.

    """

    if data.obsm["spatial"] is None:
        raise Exception("No spatial information included in AnnData object.")

    packages = ["SpaNorm", "SpatialExperiment"]
    names_to_install = [x for x in packages if not rp.isinstalled(x)]
    if len(names_to_install) > 0:
        utils.install_packages(rv.StrVector(names_to_install))
    for x in packages:
        rp.importr(x)

    params = locals()
    h.handle_r_params(params)
    r_code = f"""function(x, coord, sample.p, gene.model, adj.method, scale.factor, df.tps, lambda.a, batch, tol, 
    step.factor, maxit.nb, maxit.psi, maxn.psi, overwrite, verbose) {{
        spe <- SpatialExperiment(assays = list(counts = as.matrix(x)), spatialCoords = as.matrix(coord))
        as.data.frame(as.matrix(SpaNorm(spe,
                                        sample.p = sample.p,
                                        gene.model = gene.model,
                                        adj.method = adj.method,
                                        scale.factor = scale.factor,
                                        df.tps = df.tps,
                                        lambda.a = lambda.a,
                                        batch = batch,
                                        tol = tol,
                                        step.factor = step.factor,
                                        maxit.nb = maxit.nb,
                                        maxit.psi = maxit.psi,
                                        maxn.psi = maxn.psi,
                                        overwrite = overwrite,
                                        verbose = verbose)@assays@data@listData[["logcounts"]]))
    }}"""
    spnm = ro.r(r_code)
    
    counts = data.to_df().T
    with (ro.default_converter + rpd.converter + rnp.converter).context():
        trans = spnm(counts, data.obsm["spatial"], params["sample_p"], params["gene_model"], params["adj_method"],
                     params["scale_factor"], params["df_tps"], params["lambda_a"], params["batch"], params["tol"],
                     params["step_factor"], params["maxit_nb"], params["maxit_psi"], params["maxn_psi"],
                     params["overwrite"], params["verbose"])

    data.raw = data.copy()

    if inplace:
        data.X = trans.T
    else:
        copy = data.copy()
        copy.X = trans.T
        return copy
