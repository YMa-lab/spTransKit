import pandas as pd
import random
import scipy.stats
import numpy as np
import scanpy as sc

def get_runtime(results_dir : str, names : list) -> pd.Series:
    copy = names.copy()
    copy.remove("y")

    return pd.read_csv(results_dir + "/runtime.csv", header = 0, index_col = 0).iloc[:, 0].loc[copy]

def get_memory(results_dir : str, names : list) -> pd.Series:
    copy = names.copy()
    copy.remove("y")

    return pd.read_csv(results_dir + "/memory.csv", header = 0, index_col = 0).iloc[:, 0].loc[copy]

def create_subpanel(matrix : pd.DataFrame, genes : dict) -> list:
    subpanel = random.sample(matrix.columns.to_list(), int(len(matrix.columns.to_list()) / 4))
    for gene in genes.keys():
        if gene not in subpanel:
            subpanel.pop(0)
            subpanel.append(gene)

    return subpanel

def transform_again(counts : pd.DataFrame, names : list, num_hvg : int, subpanel : list, coordinates_dir : str) -> list:
    sub_counts = counts.loc[:, subpanel]
    sub_counts = sub_counts.loc[(sub_counts != 0).any(axis = 1), :]
    sub_counts = sub_counts.loc[(sub_counts != 0.0).any(axis = 1), :]
    
    t = []

    for name in names:
        if name == "y":
            from functions import raw
            t.append(raw.Raw(sub_counts, num_hvg))
        elif name == "y/s":
            from functions import rawsize
            t.append(rawsize.RawSize(sub_counts, num_hvg))
        elif name == "CPM":
            from functions import cpm
            t.append(cpm.CPM(sub_counts, num_hvg))
        elif name == "log(y/s + 1)":
            from functions import shiftedlog
            t.append(shiftedlog.ShiftedLog(sub_counts, num_hvg))
        elif name == "acosh(2αy/s + 1)":
            from functions import acosh
            t.append(acosh.Acosh(sub_counts, num_hvg))
        elif name == "log(y/s + 1/(4α))":
            from functions import pseudoshiftedlog
            t.append(pseudoshiftedlog.PseudoShiftedLog(sub_counts, num_hvg))
        elif name == "log(CPM + 1)":
            from functions import cpmshiftedlog
            t.append(cpmshiftedlog.CPMShiftedLog(sub_counts, num_hvg))
        elif name == "log(y/s + 1)/u":
            from functions import shiftedlogsize
            t.append(shiftedlogsize.ShiftedLogSize(sub_counts, num_hvg))
        elif name == "Analytic Pearson (no clip)":
            from functions import analyticpearson_noclip
            t.append(analyticpearson_noclip.PearsonResidual(sub_counts, num_hvg))
        elif name == "Analytic Pearson (clip)":
            from functions import analyticpearson_clip
            t.append(analyticpearson_clip.PearsonResidual(sub_counts, num_hvg))
        elif name == "scanpy Zheng":
            from functions import scanpy_zheng
            t.append(scanpy_zheng.ScanpyZheng(sub_counts, num_hvg))
        elif name == "scanpy Weinreb":
            from functions import scanpy_weinreb
            t.append(scanpy_weinreb.ScanpyWeinreb(sub_counts, num_hvg))
        elif name == "scanpy Pearson Residual":
            from functions import scanpy_pearson_residual
            t.append(scanpy_pearson_residual.ScanpyPearsonResidual(sub_counts, num_hvg))
        elif name == "DESeq2 (log)":
            from functions import deseq2_log
            t.append(deseq2_log.DESeq2(sub_counts, num_hvg))
        elif name == "DESeq2 (log1p)":
            from functions import deseq2_log1p
            t.append(deseq2_log1p.DESeq2(sub_counts, num_hvg))
        elif name == "TMM":
            from functions import tmm
            t.append(tmm.TMM(sub_counts, num_hvg))
        elif name == "Normalisr":
            from functions import normalisr
            t.append(normalisr.Normalisr(sub_counts, num_hvg))
        elif name == "PsiNorm":
            from functions import psinorm
            t.append(psinorm.PsiNorm(sub_counts, num_hvg))
        elif name == "SCTransform":
            from functions import raw
            try:
                t.append(raw.Raw(sc.read_h5ad(coordinates_dir + "/r_subcounts/sctransform.h5ad").to_df(), num_hvg, "SCTransform"))
            except(FileNotFoundError):
                t.append(None)
        elif name == "Dino":
            from functions import raw
            try:
                t.append(raw.Raw(sc.read_h5ad(coordinates_dir + "/r_subcounts/dino.h5ad").to_df(), num_hvg, "Dino"))
            except(FileNotFoundError):
                t.append(None)
        elif name == "SpaNorm":
            from functions import raw
            try:
                t.append(raw.Raw(sc.read_h5ad(coordinates_dir + "/r_subcounts/spanorm.h5ad").to_df(), num_hvg, "SpaNorm"))
            except(FileNotFoundError):
                t.append(None)
    
    return t

def sf_correspondence(trans : list[pd.DataFrame], t_two : list, names : list, results_dir : str) -> pd.Series:
    # Find size factor correspondence
    p_corr = pd.Series(None, index = names, name = "PCC")
    for t in range(0, len(trans)):
        if trans[t] is None:
            p_corr.at[names[t]] = None
        elif t_two[t] is None:
            p_corr.at[names[t]] = None
        else:
            y = t_two[t].matrix
            x = trans[t].loc[y.index.to_list(), :]
            x = (x.sum(1) / (x.sum().sum() / len(x.index.to_list()))).to_list()
            x = scipy.stats.zscore(x)
            y = (y.sum(1) / (y.sum().sum() / len(y.index.to_list()))).to_list()
            y = scipy.stats.zscore(y)
            mean_x = sum(x) / len(x)
            mean_y = sum(y) / len(y)
            x_diff = [(z - mean_x) for z in x]
            y_diff = [(z - mean_y) for z in y]
            numerator = sum(np.multiply(x_diff, y_diff).tolist())
            x_diff_sqaured = [z ** 2 for z in x_diff]
            y_diff_sqaured = [z ** 2 for z in y_diff]
            denominator = (sum(x_diff_sqaured) * sum(y_diff_sqaured)) ** 0.5
            p_corr.at[names[t]] = (numerator / denominator)

    # Save results
    p_corr.to_csv(results_dir + "/sf_correspondence.csv")
        
    # Return results
    return p_corr

def count_correspondence(trans : list[pd.DataFrame], t_two : list, names : list, genes : dict, results_dir : str) -> pd.Series:
    # Find transformed count correspondence
    gene_list = list(genes.keys())
    p_corr = pd.DataFrame(None, index = names, columns = gene_list)
    for gene in gene_list:
        for t in range(0, len(trans)):
            if trans[t] is None:
                p_corr.at[names[t], gene] = None
            elif t_two[t] is None:
                p_corr.at[names[t], gene] = None
            else:
                y = t_two[t].matrix
                x = trans[t].loc[y.index.to_list(), gene].to_list()
                y = t_two[t].matrix.loc[:, gene].to_list()
                mean_x = sum(x) / len(x)
                mean_y = sum(y) / len(y)
                x_diff = [(z - mean_x) for z in x]
                y_diff = [(z - mean_y) for z in y]
                numerator = sum(np.multiply(x_diff, y_diff).tolist())
                x_diff_sqaured = [z ** 2 for z in x_diff]
                y_diff_sqaured = [z ** 2 for z in y_diff]
                denominator = (sum(x_diff_sqaured) * sum(y_diff_sqaured)) ** 0.5
                if denominator == 0.0:
                    p_corr.at[names[t], gene] = None
                else:
                    p_corr.at[names[t], gene] = numerator / denominator

    # Save results
    p_corr.to_csv(results_dir + "/count_correspondence.csv")
    
    # Return results
    return p_corr.mean(1)
