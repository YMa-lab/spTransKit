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

import helpers as h


# Size Factor-Based Transformations

def size(matrix : sc.AnnData, inplace : bool = False):
    """ This function applies the y/s transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False

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
    matrix
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    if inplace:
        matrix.X = np.divide(matrix.X, h.cpm(matrix.X))
    else:
        copy = matrix.copy()
        copy.X = np.divide(matrix.X, h.cpm(matrix.X))
        return copy
    

def zheng(matrix : sc.AnnData, inplace : bool = False):
    """ This function applies the scanpy transformation to the gene count matrix according to Zheng et al., 2016.
    
    Parameters
    ----------
    matrix
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    if inplace:
        sc.pp.normalize_total(matrix, key_added = "n_counts_all")
        sc.pp.log1p(matrix.X)
        sc.pp.scale(matrix.X)
    else:
        copy = matrix.copy()
        sc.pp.normalize_total(copy, key_added = "n_counts_all")
        sc.pp.log1p(copy.X)
        sc.pp.scale(copy.X)
        return copy
    

def tmm(matrix : sc.AnnData,
        logratio_trim : float = 0.3,
        sum_trim : float = 0.05,
        inplace : bool = False):
    """ This function applies the TMM transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
    logratio_trim
        float, fraction of data to trim from the ends of ranked log ratios, default = 0.3
    sum_trim
        float, fraction of data to trim from the ends of ranked sums, default = 0.05
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    np.seterr(all = "ignore")

    m = pd.DataFrame(matrix.X)

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
    matrix
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """
    np.seterr(all = "ignore")

    m = pd.DataFrame(matrix.X)

    log_counts = pd.DataFrame(np.log(m), index = m.index, columns = m.columns)
    log_counts.replace(float("-inf"), np.nan, inplace = True)
    log_means = log_counts.mean(0, skipna = True)
    filtered_genes = log_means.loc[log_means > 0].index.to_list()
    log_ratios = log_counts.loc[:, filtered_genes].sub(log_means.loc[filtered_genes], axis = 1)

    log_medians = log_ratios.median(axis = 1)
    log_medians.fillna(0, inplace = True)
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
    matrix
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """
    
    shift = np.log1p((matrix.X / h.size_factor(matrix.X)))

    if inplace:
        matrix.X = shift
    else:
        copy = matrix.copy()
        copy.X = shift
        return copy
    

def cpm_shifted_log(matrix : sc.AnnData, inplace : bool = False):
    """ This function applies the log(CPM + 1) transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    shift = np.log1p(matrix.X / h.cpm(matrix.X))
    
    if inplace:
        matrix.X = shift
    else:
        copy = matrix.copy()
        copy.X = shift
        return copy


def shifted_log_size(matrix : sc.AnnData, inplace : bool = False):
    """ This function applies the log(y/s + 1)/u transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    shift = np.log1p(matrix.X / h.size_factor(matrix.X))

    if inplace:
        matrix.X = shift / h.size_factor(shift)
    else:
        copy = matrix.copy()
        copy.X = shift / h.size_factor(shift)
        return copy


def acosh(matrix : sc.AnnData, alpha : float = 0.05, inplace : bool = False):
    """ This function applies the acosh(2αy/s + 1) transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
    alpha
        float, overdispersion hyperparamter, default = 0.05
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    size = matrix.X / h.size_factor(matrix.X)
    shift = (1 / np.sqrt(alpha)) * np.arccosh((size * alpha * 2) + 1)

    if inplace:
        matrix.X = shift
    else:
        copy = matrix.copy()
        copy.X = shift
        return copy
    

def pseudo_shifted_log(matrix : sc.AnnData, alpha : float = 0.05, inplace : bool = False):
    """ This function applies the log(y/s + 1/4α) transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
    alpha
        float, overdispersion hyperparamter, default = 0.05
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    shift = (1 / np.sqrt(alpha)) * np.log1p(4 * alpha * (matrix.X / h.size_factor(matrix.X)))

    if inplace:
        matrix.X = shift
    else:
        copy = matrix.copy()
        copy.X = shift
        return copy
    

# Model-Based Transformations

def analytic_pearson_noclip(matrix : sc.AnnData, alpha : float = 0.05, inplace : bool = False):
    """ This function applies the Analytic Pearson (no clip) transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
    alpha
        float, overdispersion hyperparamter, default = 0.05
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False

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


def analytic_pearson_clip(matrix : sc.AnnData, alpha : float = 0.05, inplace : bool = False):
    """ This function applies the Analytic Pearson (clip) transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
    alpha
        float, overdispersion hyperparamter, default = 0.05
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    nb_matrix = np.zeros(matrix.X.shape)
    nb_matrix = np.add(nb_matrix, h.n(matrix.X))
    nb_matrix = np.multiply(nb_matrix, h.p(matrix.X))

    var_matrix = nb_matrix.copy()
    var_matrix = np.add(((var_matrix ** 2) * alpha), nb_matrix) ** 0.5

    final_matrix = np.divide((np.subtract(matrix.X, nb_matrix)), var_matrix)
    
    clip = matrix.X.shape[0] ** 0.5
    final_matrix = np.where(final_matrix > clip, clip, final_matrix)
    final_matrix = np.where(final_matrix < -clip, -clip, final_matrix)
    
    if inplace:
        matrix.X = final_matrix
    else:
        copy = matrix.copy()
        copy.X = final_matrix
        return copy


def sc_pearson(matrix : sc.AnnData, 
               theta : float = 100,
               clip : float = None,
               check_values : bool = True,
               layer : str = None,
               inplace : bool = False):
    """ This function applies the scanpy Pearson residual transformation to the gene count matrix.
    
    Parameters
    ----------
    matrix
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
    theta
        float,
    clip
        float,
    check_values
        bool,
    layer
        str,
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    if inplace:
        sc.experimental.pp.normalize_pearson_residuals(matrix, theta = theta, clip = clip, check_values = check_values, 
                                                       layer = layer)
    else:
        copy = matrix.copy()
        sc.experimental.pp.normalize_pearson_residuals(copy)
        return copy


def normalisr(matrix : sc.AnnData, 
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
    matrix
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
    normalize,
        bool,
    nth
        int,
    ntot
        int,
    varscale
        float,
    seed
        int,
    lowmem
        bool,
    nocov
        bool,
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    norm = n.lcpm(np.transpose(matrix.X), normalize = normalize, nth = nth, ntot = ntot, varscale = varscale, seed = seed,
                  lowmem = lowmem, nocov = nocov)[0]

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
    matrix
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    n = matrix.shape[0]
    m = (matrix.X + 1).min(1)
    log_counts = np.log(matrix.X + 1)
    log_m = np.log(m)
    log_m = np.reshape(log_m, (matrix.shape[0], 1))
    sf = np.subtract(log_counts, log_m).sum(1) / n
    sf = np.reshape(sf, (matrix.shape[0], 1))
    
    if inplace:
        matrix.X = np.divide(matrix.X, sf)
    else:
        copy = matrix.copy()
        copy.X = np.divide(matrix.X, sf)
        return copy


def sctransform(matrix : sc.AnnData, 
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
                min_variance : float | Literal["umi_median", "model_median", "model_mean"] = "-Inf",
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
    matrix
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
    cell_attr
        pd.DataFrame,
    latent_var
        list[str],
    batch_var
        list[str],
    latent_var_nonreg
        list[str],
    n_genes
        int,
    n_cells
        int,
    method
        str,
    do_regularize
        bool,
    theta_regularization
        str,
    res_clip_range
        list[str],
    bin_size
        int,
    min_cells
        int,
    residual_type
        str,
    return_cell_attr
        bool,
    return_gene_attr
        bool,
    return_corrected_umi
        bool,
    min_variance
        float or str,
    bw_adjust
        float,
    gmean_eps
        float,
    theta_estimation_fun
        str,
    theta_given
        list or float,
    exclude_poisson
        bool,
    use_geometric_mean
        bool,
    use_geometric_mean_offset
        bool,
    fix_intercept
        bool,
    fix_slope
        bool,
    scale_factor
        float,
    vst_flavor
        str,
    verbosity
        int,
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

    packages = ["sctransform"]
    names_to_install = [x for x in packages if not rp.isinstalled(x)]
    if len(names_to_install) > 0:
        utils.install_packages(rv.StrVector(names_to_install))
    for x in packages:
        rp.importr(x)

    params = locals()
    h.handle_r_params(params)
    if res_clip_range is None:
        params["res_clip_range"] = ro.FloatVector([-(matrix.shape[0] ** 0.5), matrix.shape[0] ** 0.5])
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

    counts = matrix.to_df().T
    with (ro.default_converter + rpd.converter).context():
        trans = vst(counts, params["cell_attr"], params["latent_var"], params["batch_var"], params["latent_var_nonreg"],
                    params["n_genes"], params["n_cells"], params["method"], params["do_regularize"],
                    params["theta_regularization"], params["res_clip_range"], params["bin_size"], params["min_cells"],
                    params["residual_type"], params["return_cell_attr"], params["return_gene_attr"],
                    params["return_corrected_umi"], params["min_variance"], params["bw_adjust"], params["gmean_eps"],
                    params["theta_estimation_fun"], params["theta_given"], params["exclude_poisson"],
                    params["use_geometric_mean"], params["use_geometric_mean_offset"], params["fix_intercept"],
                    params["fix_slope"], params["scale_factor"], params["vst_flavor"], params["verbosity"])
                
    if inplace:
        matrix = sc.AnnData(trans.T)
    else:
        return sc.AnnData(trans.T)


def dino(matrix : sc.AnnData, 
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
    matrix
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
    nCores
        int,
    prec
        int,
    minNZ
        int,
    nSubGene
        int,
    nSubCell
        int,
    depth
        list[float],
    slope
        float,
    minSlope
        float,
    maxSlope
        float,
    clusterSlope
        bool,
    returnMeta
        bool,
    doRQS
        bool,
    emPar
        dict,
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

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

    counts = matrix.to_df().T
    with (ro.default_converter + rpd.converter).context():
        trans = dino(counts, params["nCores"], params["prec"], params["minNZ"], params["nSubGene"], params["nSubCell"], 
                     params["depth"], params["slope"], params["minSlope"], params["maxSlope"], params["clusterSlope"],
                     params["returnMeta"], params["doRQS"], params["emPar"])

    if inplace:
        matrix.X = trans.T
    else:
        copy = matrix.copy()
        copy.X = trans.T
        return copy


# Spatially Aware Transformations

def spanorm(matrix : sc.AnnData,
            coordinates : pd.DataFrame,
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
    matrix
        sc.AnnData, N x G gene count matrix. Rows correspond to spatial locations and columns correspond to genes.
    coordinates
        pd.DataFrame, N x 2 coordinate matrix. Rows correspond to spatial locations and columns correspond to the x and
        y spatial coordinates.
    sample_p
        float,
    gene_model
        str,
    adj_method
        str,
    scale_factor
        float,
    df_tps
        float,
    lambda_a
        float,
    batch
        pd.DataFrame,
    tol
        float,
    step_factor
        float,
    maxit_nb
        int,
    maxit_psi
        int,
    maxn_psi
        int,
    overwrite
        bool,
    verbose
        bool,
    inplace
        bool, whether to transform the matrix within the original AnnData object, default = False

    Returns
    ----------
    copy
        sc.AnnData, copy of the transformed gene count matrix if inplace = False

    """

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
    
    counts = matrix.to_df().T
    with (ro.default_converter + rpd.converter + rnp.converter).context():
        trans = spnm(counts, coordinates, params["sample_p"], params["gene_model"], params["adj_method"],
                     params["scale_factor"], params["df_tps"], params["lambda_a"], params["batch"], params["tol"],
                     params["step_factor"], params["maxit_nb"], params["maxit_psi"], params["maxn_psi"],
                     params["overwrite"], params["verbose"])

    if inplace:
        matrix.X = trans.T
    else:
        copy = matrix.copy()
        copy.X = trans.T
        return copy
