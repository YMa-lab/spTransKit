import scanpy
import pandas
import random

import transformations

def get_feat(t : list, num_feat : int):
    x = []
    y = []
    z = []
    for i in range(1, len(t)):
        x.append(t[i].name)
        y.append(t[i].runtime)
        z.append(t[i].memory)
    return pandas.Series(y, index = x, name = num_feat), pandas.Series(z, index = x, name = num_feat)

def sub_spot(matrix : pandas.DataFrame, num_spot : int) -> pandas.DataFrame:
    sub = matrix.loc[random.sample(matrix.index.to_list(), num_spot), :]
    sub = sub.loc[(sub != 0).any(axis = 1), :]
    sub = sub.loc[:, (sub != 0).any(axis = 0)]
    return sub

def compeff_spot():
    matrix = scanpy.read_h5ad("/oscar/data/yma16/Project/spTransform/0.NormalizedCounts/10x_Xenium_Prime/Human_Prostate/norm_counts/raw.h5ad").to_df()
    names = ["y", "y/s", "CPM", "log(y/s + 1)", "acosh(2αy/s + 1)", "log(y/s + 1/(4α))", "log(CPM + 1)", "log(y/s + 1)/u", "Analytic Pearson (no clip)", "Analytic Pearson (clip)", "scanpy Zheng", "scanpy Pearson Residual", "DESeq2 (log)", "DESeq2 (log1p)", "TMM", "Normalisr", "PsiNorm"]
    sub_1000_r, sub_1000_m = get_feat(transformations.transform(sub_spot(matrix, 1000), names, 10), 1000)
    sub_5000_r, sub_5000_m = get_feat(transformations.transform(sub_spot(matrix, 5000), names, 10), 5000)
    sub_10000_r, sub_10000_m = get_feat(transformations.transform(sub_spot(matrix, 10000), names, 10), 10000)
    sub_25000_r, sub_25000_m = get_feat(transformations.transform(sub_spot(matrix, 25000), names, 10), 25000)
    sub_50000_r, sub_50000_m = get_feat(transformations.transform(sub_spot(matrix, 50000), names, 10), 50000)
    sub_100000_r, sub_100000_m = get_feat(transformations.transform(sub_spot(matrix, 100000), names, 10), 100000)
    sub_125000_r, sub_125000_m = get_feat(transformations.transform(sub_spot(matrix, 125000), names, 10), 125000)
    sub_150000_r, sub_150000_m = get_feat(transformations.transform(sub_spot(matrix, 150000), names, 10), 150000)
    sub_175000_r, sub_175000_m = get_feat(transformations.transform(sub_spot(matrix, 175000), names, 10), 175000)
    sub_max_r, sub_max_m = get_feat(transformations.transform(matrix, names, 10), len(matrix.index.to_list()))
    pandas.concat([sub_1000_r, sub_5000_r, sub_10000_r, sub_25000_r, sub_50000_r, sub_100000_r, sub_125000_r, sub_150000_r, sub_175000_r, sub_max_r], axis = 1).to_csv("/oscar/data/yma16/Project/spTransform/1.Evaluations/10x_Xenium_Prime/Human_Prostate/july_results/spot_runtime_py.csv")
    pandas.concat([sub_1000_m, sub_5000_m, sub_10000_m, sub_25000_m, sub_50000_m, sub_100000_m, sub_125000_m, sub_150000_m, sub_175000_m, sub_max_m], axis = 1).to_csv("/oscar/data/yma16/Project/spTransform/1.Evaluations/10x_Xenium_Prime/Human_Prostate/july_results/spot_memory_py.csv")

def sub_gene(matrix : pandas.DataFrame, num_gene : int) -> pandas.DataFrame:
    sub = matrix.loc[:, random.sample(matrix.columns.to_list(), num_gene)]
    sub = sub.loc[(sub != 0).any(axis = 1), :]
    sub = sub.loc[:, (sub != 0).any(axis = 0)]
    return sub

def compeff_gene():
    matrix = scanpy.read_h5ad("/oscar/data/yma16/Project/spTransform/0.NormalizedCounts/Open-ST/Human_MLN_3/norm_counts/raw.h5ad").to_df()
    names = ["y", "y/s", "CPM", "log(y/s + 1)", "acosh(2αy/s + 1)", "log(y/s + 1/(4α))", "log(CPM + 1)", "log(y/s + 1)/u", "Analytic Pearson (no clip)", "Analytic Pearson (clip)", "scanpy Zheng", "scanpy Pearson Residual", "DESeq2 (log)", "DESeq2 (log1p)", "TMM", "Normalisr", "PsiNorm"]
    sub_1000_r, sub_1000_m = get_feat(transformations.transform(sub_gene(matrix, 1000), names, 10), 1000)
    sub_5000_r, sub_5000_m = get_feat(transformations.transform(sub_gene(matrix, 5000), names, 10), 5000)
    sub_10000_r, sub_10000_m = get_feat(transformations.transform(sub_gene(matrix, 10000), names, 10), 10000)
    sub_15000_r, sub_15000_m = get_feat(transformations.transform(sub_gene(matrix, 15000), names, 10), 15000)
    sub_20000_r, sub_20000_m = get_feat(transformations.transform(sub_gene(matrix, 20000), names, 10), 20000)
    sub_25000_r, sub_25000_m = get_feat(transformations.transform(sub_gene(matrix, 25000), names, 10), 25000)
    sub_30000_r, sub_30000_m = get_feat(transformations.transform(sub_gene(matrix, 30000), names, 10), 30000)
    sub_35000_r, sub_35000_m = get_feat(transformations.transform(sub_gene(matrix, 35000), names, 10), 35000)
    sub_max_r, sub_max_m = get_feat(transformations.transform(matrix, names, 10), len(matrix.columns.to_list()))
    pandas.concat([sub_1000_r, sub_5000_r, sub_10000_r, sub_15000_r, sub_20000_r, sub_25000_r, sub_30000_r, sub_35000_r, sub_max_r], axis = 1).to_csv("/oscar/data/yma16/Project/spTransform/1.Evaluations/Open-ST/Human_MLN_3/july_results/gene_runtime_py.csv")
    pandas.concat([sub_1000_m, sub_5000_m, sub_10000_m, sub_15000_m, sub_20000_m, sub_25000_m, sub_30000_m, sub_35000_m, sub_max_m], axis = 1).to_csv("/oscar/data/yma16/Project/spTransform/1.Evaluations/Open-ST/Human_MLN_3/july_results/gene_memory_py.csv")

if __name__ == "__main__":
    compeff_spot()
    compeff_gene()