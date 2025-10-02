import argparse
import os
import pandas as pd
import scanpy as sc

def transform(counts : pd.DataFrame, names : list, num_hvg : int) -> list:
    t = []

    for name in names:
        if name == "y":
            from functions import raw
            t.append(raw.Raw(counts, num_hvg))
        elif name == "y/s":
            from functions import rawsize
            t.append(rawsize.RawSize(counts, num_hvg))
        elif name == "CPM":
            from functions import cpm
            t.append(cpm.CPM(counts, num_hvg))
        elif name == "log(y/s + 1)":
            from functions import shiftedlog
            t.append(shiftedlog.ShiftedLog(counts, num_hvg))
        elif name == "acosh(2αy/s + 1)":
            from functions import acosh
            t.append(acosh.Acosh(counts, num_hvg))
        elif name == "log(y/s + 1/(4α))":
            from functions import pseudoshiftedlog
            t.append(pseudoshiftedlog.PseudoShiftedLog(counts, num_hvg))
        elif name == "log(CPM + 1)":
            from functions import cpmshiftedlog
            t.append(cpmshiftedlog.CPMShiftedLog(counts, num_hvg))
        elif name == "log(y/s + 1)/u":
            from functions import shiftedlogsize
            t.append(shiftedlogsize.ShiftedLogSize(counts, num_hvg))
        elif name == "Analytic Pearson (no clip)":
            from functions import analyticpearson_noclip
            t.append(analyticpearson_noclip.PearsonResidual(counts, num_hvg))
        elif name == "Analytic Pearson (clip)":
            from functions import analyticpearson_clip
            t.append(analyticpearson_clip.PearsonResidual(counts, num_hvg))
        elif name == "scanpy Zheng":
            from functions import scanpy_zheng
            t.append(scanpy_zheng.ScanpyZheng(counts, num_hvg))
        elif name == "scanpy Weinreb":
            from functions import scanpy_weinreb
            t.append(scanpy_weinreb.ScanpyWeinreb(counts, num_hvg))
        elif name == "scanpy Pearson Residual":
            from functions import scanpy_pearson_residual
            t.append(scanpy_pearson_residual.ScanpyPearsonResidual(counts, num_hvg))
        elif name == "DESeq2 (log)":
            from functions import deseq2_log
            t.append(deseq2_log.DESeq2(counts, num_hvg))
        elif name == "DESeq2 (log1p)":
            from functions import deseq2_log1p
            t.append(deseq2_log1p.DESeq2(counts, num_hvg))
        elif name == "TMM":
            from functions import tmm
            t.append(tmm.TMM(counts, num_hvg))
        elif name == "Normalisr":
            from functions import normalisr
            t.append(normalisr.Normalisr(counts, num_hvg))
        elif name == "PsiNorm":
            from functions import psinorm
            t.append(psinorm.PsiNorm(counts, num_hvg))

    return t

def save_counts(t_one : list, files : list, dir : str):
    for i in range(0, len(t_one)):
        sc.AnnData(t_one[i].matrix).write_h5ad(dir + "/" + files[i] + ".h5ad", compression = "gzip")

def save_hvg(t_one : list, files : list, hvg_dir : str):
    for i in range(0, len(t_one)):
        pd.Series(t_one[i].hvgenes, name = "HVG").to_csv(hvg_dir + "/" + files[i] + ".csv")

def save_runtime(t_one : list, names : list, results_dir : str):
    runtime = pd.Series(0.0, name = "Runtime (s)", index = names)
    for t in t_one:
        runtime[t.name] = t.runtime

    runtime.to_csv(results_dir + "/runtime_p.csv")

def save_memory(t_one : list, names : list, results_dir : str):
    memory = pd.Series(0.0, name = "Memory (MB)", index = names)
    for t in t_one:
        memory[t.name] = t.memory

    memory.to_csv(results_dir + "/memory_p.csv")

if __name__ == "__main__":
    # Parse user input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_id", type = str, help = "sample identifier")
    parser.add_argument("--technology", type = str, help = "type of ST technology")
    parser.add_argument("--num_hvg", type = int, help = "number of highly variable genes to be selected")
    params = parser.parse_args()

    # Create directories
    dir = "/oscar/data/yma16/Project/spTransform/"

    coordinates_dir = dir + "0.NormalizedCounts/" + params.technology + "/" + params.sample_id
    if not os.path.exists(coordinates_dir):
        os.makedirs(coordinates_dir)

    norm_counts_dir = coordinates_dir + "/norm_counts"
    if not os.path.exists(norm_counts_dir):
        os.makedirs(norm_counts_dir)

    hvg_dir = dir + "1.Evaluations/" + params.technology + "/" + params.sample_id + "/july_hvg"
    if not os.path.exists(hvg_dir):
        os.makedirs(hvg_dir)

    results_dir = dir + "1.Evaluations/" + params.technology + "/" + params.sample_id + "/july_results"
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    # Setup data
    names = ["y", "y/s", "CPM", "log(y/s + 1)", "acosh(2αy/s + 1)", "log(y/s + 1/(4α))", "log(CPM + 1)", "log(y/s + 1)/u", "Analytic Pearson (no clip)", "Analytic Pearson (clip)", "scanpy Zheng", "scanpy Pearson Residual", "DESeq2 (log)", "DESeq2 (log1p)", "TMM", "Normalisr", "PsiNorm"]
    files = ["raw", "raw_size", "cpm", "shifted_log", "acosh", "pseudo_shifted_log", "cpm_shifted_log", "shifted_log_size", "analytic_pearson_residual_noclip", "analytic_pearson_residual_clip", "scanpy_zheng", "scanpy_pearson_residual", "deseq2_log", "deseq2_log1p", "tmm", "normalisr", "psinorm"]
    counts = sc.read_h5ad(norm_counts_dir + "/raw.h5ad").to_df()
    
    # 2. Transformation
    t_one = transform(counts, names, params.num_hvg)
    save_counts(t_one, files, norm_counts_dir)
    save_hvg(t_one, files, hvg_dir)
    save_runtime(t_one, names, results_dir)
    save_memory(t_one, names, results_dir)