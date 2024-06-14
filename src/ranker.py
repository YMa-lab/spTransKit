import scipy.stats
import pandas
import numpy

def rank_runtime(transformations : list) -> list:
    y = []
    for t in transformations:
        y.append(t.runtime)
    max_rt = max(y)
    return [1 - (z / max_rt) for z in y]

def rank_memory(transformations : list) -> list:
    y = []
    for t in transformations:
        y.append(t.memory)
    max_mem = max(y)
    return [1 - (z / max_mem) for z in y]

def rank_gene_mean_var(transformations : list) -> list:
    ks_values = []
    for t in transformations:
        ks_values.append(scipy.stats.kstest(scipy.stats.zscore(t.matrix.var(0)), "norm")[1])
    return ks_values

def morans_i(weight_matrix : pandas.DataFrame, gene_dist : pandas.Series) -> float:
    n = gene_dist.size
    mean = gene_dist.mean()
    i_dist = gene_dist - mean
    i_dist_squared = i_dist ** 2
    j_dist = gene_dist - mean
    numerator = pandas.DataFrame(weight_matrix, copy = True)
    numerator = numerator.multiply(i_dist, axis = 1)
    numerator = numerator.multiply(j_dist, axis = 0)
    weight_sum = weight_matrix.sum(0).sum()
    if weight_sum * i_dist_squared.sum() == 0.0:
        moran = 0.0
    else:
        moran = (n * numerator.sum(0).sum()) / (weight_sum * i_dist_squared.sum())
    return moran

def rank_marker_genes(transformations : list, name : list, genes : dict, weight_matrix : pandas.DataFrame, results_dir : str) -> list:
    gene_dist = pandas.DataFrame(0.0, index = name, columns = genes.keys())
    for gene in genes.keys():
        if genes[gene] == "normal":
            for t in transformations:
                gene_dist.loc[t.name, gene] = scipy.stats.kstest(scipy.stats.zscore(t.matrix.loc[:, gene]), "norm")[1]
        elif genes[gene] == "spatial":
            for t in transformations:
                gene_dist.loc[t.name, gene] = morans_i(weight_matrix, t.matrix.loc[:, gene])

    gene_dist.to_csv(results_dir + "/marker_gene_distributions.csv")
    
    return gene_dist.mean(1).to_list()

def rank_hvg_similarity(transformations : list) -> list:
    y = []
    raw = transformations[0].hvgenes
    for t in transformations:
        trans = t.hvgenes
        intersection = len(set(raw).intersection(set(trans)))
        union = len(set(raw).union(set(trans)))
        jaccard = intersection / union
        y.append(jaccard)
    return y

def rank_pathway_enrichment(enrichment : pandas.DataFrame) -> list:
    enr = enrichment.sum(0)
    max_enr = max(enr)
    if max_enr == 0.0:
        return [0.0 for z in enr]
    else:
        return [(z / max_enr) for z in enr]
    
def rank_sf_correspondence(t_one : list, t_two : list, name : list, results_dir : str) -> list:
    p_corr = []
    for t in range(0, len(t_one)):
        y = t_two[t].matrix
        x = t_one[t].matrix.loc[y.index.to_list(), :]
        x = (x.sum(1) / (x.sum().sum() / len(x.index.to_list()))).to_list()
        x = scipy.stats.zscore(x)
        y = (y.sum(1) / (y.sum().sum() / len(y.index.to_list()))).to_list()
        y = scipy.stats.zscore(y)
        mean_x = sum(x) / len(x)
        mean_y = sum(y) / len(y)
        x_diff = [(z - mean_x) for z in x]
        y_diff = [(z - mean_y) for z in y]
        numerator = sum(numpy.multiply(x_diff, y_diff).tolist())
        x_diff_sqaured = [z ** 2 for z in x_diff]
        y_diff_sqaured = [z ** 2 for z in y_diff]
        denominator = (sum(x_diff_sqaured) * sum(y_diff_sqaured)) ** 0.5
        p_corr.append(numerator / denominator)
    
    pandas.Series(p_corr, index = name, name = "Pearson Correlation").to_csv(results_dir + "/sf_correspondence.csv")

    min_corr = min(p_corr)
    if min_corr < 0:
        p_corr = [(z - min_corr) for z in p_corr]
    max_corr = max(p_corr)
    p_corr = [(z / max_corr) for z in p_corr]
    return p_corr

def rank_count_correspondence(t_one : list, t_two : list, genes : dict, name : list, results_dir : str) -> list:
    p_corr = []
    gene_list = list(genes.keys())
    for i in range(0, len(gene_list)):
        p_corr.append([])
        for t in range(0, len(t_one)):
            y = t_two[t].matrix
            x = t_one[t].matrix.loc[y.index.to_list(), gene_list[i]].to_list()
            y = t_two[t].matrix.loc[:, gene_list[i]].to_list()
            mean_x = sum(x) / len(x)
            mean_y = sum(y) / len(y)
            x_diff = [(z - mean_x) for z in x]
            y_diff = [(z - mean_y) for z in y]
            numerator = sum(numpy.multiply(x_diff, y_diff).tolist())
            x_diff_sqaured = [z ** 2 for z in x_diff]
            y_diff_sqaured = [z ** 2 for z in y_diff]
            denominator = (sum(x_diff_sqaured) * sum(y_diff_sqaured)) ** 0.5
            p_corr[i].append(numerator / denominator)

        pandas.Series(p_corr[i], index = name, name = "Pearson Correlation").to_csv(results_dir + "/" + gene_list[i] + "_count_correspondence.csv")

        min_corr = min(p_corr[i])
        if min_corr < 0:
            p_corr[i] = [(z - min_corr) for z in p_corr[i]]
        max_corr = max(p_corr[i])
        p_corr[i] = [(z / max_corr) for z in p_corr[i]]
    pcavg = pandas.DataFrame(p_corr).mean(0).to_list()
    return pcavg

def rank_corr_matrix(transformations : list, name : list, results_dir : str) -> list:
    corr_matrix_sims = []
    raw = transformations[0].matrix.loc[:, transformations[0].hvgenes].corr()
    n = raw.shape[0] * raw.shape[1]
    for t in range(0, len(transformations)):
       squared_err = transformations[t].matrix.loc[:, transformations[t].hvgenes].corr().sub(raw) ** 2
       sse = squared_err.sum(0).sum()
       rmse = (sse / n) ** 0.5
       corr_matrix_sims.append(rmse)

    pandas.Series(corr_matrix_sims, index = name, name = "RMSE").to_csv(results_dir + "/corr_matrix.csv")

    max_rmse = max(corr_matrix_sims)
    corr_matrix_sims = [(1 - (z / max_rmse)) for z in corr_matrix_sims]
    return corr_matrix_sims
