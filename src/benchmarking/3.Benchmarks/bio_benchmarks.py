import pandas as pd
import numpy as np
import scanpy as sc
import sklearn.neighbors
import reactome2py.analysis

import b_helpers as h

def hvg_similarity(hvgs : list[list], names : list, results_dir : str) -> pd.Series:
    # Find similarity
    sims = pd.DataFrame(0.0, index = names, columns = names)
    for i in range(0, len(hvgs)):
        i_hvg = hvgs[i]
        for j in range(0, len(hvgs)):
            j_hvg = hvgs[j]
            if i_hvg is None or j_hvg is None:
                sims.at[names[i], names[j]] = None
            else:
                intersection = len(set(i_hvg).intersection(set(j_hvg)))
                union = len(set(i_hvg).union(set(j_hvg)))
                jaccard = intersection / union
                sims.at[names[i], names[j]] = jaccard

    # Save results
    sims.to_csv(results_dir + "/hvg_similarity.csv")
    
    # Return results
    similarity = sims.iloc[:, 0]
    consistency = sims.iloc[:, 1:len(names)].mean(1)
    consistency.name = "Consistency"
    return similarity, consistency

def gene_correlation(trans : list[pd.DataFrame], hvgs : list[list], names : list, results_dir : str) -> pd.Series:
    # Find RMSE
    corr_matrix_rmse = pd.Series(0.0, index = names, name = "RMSE")
    raw = trans[0].loc[:, hvgs[0]]
    for i in range(1, len(trans)):
        if trans[i] is None:
            corr_matrix_rmse.at[names[i]] = None
        else:
            t = trans[i].loc[:, hvgs[i]]
            shared_hvg = list(set(raw.columns.to_list()).intersection(set(t.columns.to_list())))
            squared_err = t.loc[:, shared_hvg].corr().sub(raw.loc[:, shared_hvg].corr()) ** 2
            sse = squared_err.sum(0).sum()
            rmse = (sse / (len(shared_hvg) ** 2)) ** 0.5
            corr_matrix_rmse.at[names[i]] = rmse

    # Save results
    corr_matrix_rmse.to_csv(results_dir + "/corr_matrix.csv")
    
    # Return results
    return corr_matrix_rmse


def knn_similarity(trans : list[pd.DataFrame], hvgs : list[list], files : list, names : list, k : int, knn_dir : str, results_dir : str) -> pd.Series:
    # Find neighbors
    neighbors = []
    for i in range(0, len(trans)):
        if trans[i] is None:
            neighbors.append(None)
        else:
            hvg_matrix = trans[i].loc[:, hvgs[i]]
            ball_tree = sklearn.neighbors.BallTree(hvg_matrix)
            cell_to_neighbors = {}
            for cell in hvg_matrix.index.to_list():
                dist, ind = ball_tree.query(pd.DataFrame(hvg_matrix.loc[cell, :]).T, k = k + 1)
                cell_to_neighbors[cell] = hvg_matrix.iloc[ind.tolist()[0], :].index.to_list()
                if cell in cell_to_neighbors[cell]:
                    cell_to_neighbors[cell].remove(cell)
                else:
                    cell_to_neighbors[cell].pop(-1)
            # Save neighbors
            pd.DataFrame().from_dict(data = cell_to_neighbors, orient = "index").to_csv(knn_dir + "/" + files[i] + ".csv")
            neighbors.append(cell_to_neighbors)

    # Find similarity
    raw = neighbors[0]
    knn_sims = []
    for t in range(0, len(neighbors)):
        if neighbors[t] is None:
            knn_sims.append(None)
        else:
            total = 0.0
            for cell in raw.keys():
                sim = len(set(neighbors[t][cell]).intersection(set(raw[cell]))) / len(set(neighbors[t][cell]).union(set(raw[cell])))
                total = total + sim
            knn_sims.append(total / len(raw.keys()))
    sims = pd.Series(knn_sims, index = names, name = "Similarity")

    # Save results
    sims.to_csv(results_dir + "/knn_similarity.csv")
    
    # Return results
    return sims


def gene_distributions(coordinates : pd.DataFrame, genes : dict, trans : list[pd.DataFrame], names : list, results_dir : str) -> pd.Series:
    # Find Moran's I
    wm = h.weight_matrix(coordinates)
    gene_dist = pd.DataFrame(0.0, index = names, columns = genes.keys())
    for gene in genes.keys():
        for i in range(0, len(trans)):
            if trans[i] is None:
                gene_dist.at[names[i], gene] = None
            else:
                gene_dist.at[names[i], gene] = h.morans_i(wm, trans[i].loc[:, gene])

    # Save results
    gene_dist.to_csv(results_dir + "/marker_gene_distribution.csv")
    
    # Return results
    for gene in genes.keys():
        if genes[gene] == "spatial":
            gene_dist.loc[:, gene] = (gene_dist.loc[:, gene] - 1).abs()
        elif genes[gene] == "normal":
            gene_dist.loc[:, gene] = gene_dist.loc[:, gene].abs()
    gene_dist_mean = gene_dist.mean(1)
    gene_dist_mean.name = "Deviation"
    return gene_dist_mean

def pathway_enrichment(pwys : list, hvgs : list[list], names : list, species : str, results_dir : str) -> pd.Series:
    # Find marker pathways
    y = pd.DataFrame(1.0, index = pwys, columns = names)
    i = 0
    for i in range(0, len(hvgs)):
        if hvgs[i] is None:
            y[names[i]] = pd.Series(None, index = pwys)
        else:
            hvg = ""
            for gene in hvgs[i]:
                hvg = hvg + gene + ","
            result = reactome2py.analysis.identifiers(ids = hvg, species = species)
            token = result["summary"]["token"]
            marker_pathways = ""
            for pwy in pwys:
                marker_pathways = marker_pathways + pwy + ","
            pathways = reactome2py.analysis.token_pathways_result(token, marker_pathways, species = species)
            pathways_series = pd.Series(1.0, index = pwys)
            for pwy in pathways:
                pathways_series[pwy["stId"]] = float(pwy["entities"]["pValue"])
            y[names[i]] = pathways_series
            i = i + 1
    enr = (-pd.DataFrame(np.log10(y)))
    
    # Save results
    enr.to_csv(results_dir + "/pathway_enrichment.csv")

    # Return results
    enr = enr.where(enr > 0, 1)
    enr = enr.sum(0) / len(pwys)
    enr.name = "Proportion"
    return enr
