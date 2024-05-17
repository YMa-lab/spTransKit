import matplotlib.pyplot as plt
import seaborn
import pandas
import scipy.stats

import ranker

def plot_runtime(transformations : list, sample_id : str):
    x = []
    y = []
    for i in range(1, len(transformations)):
        x.append(transformations[i].name)
        y.append(transformations[i].runtime)
    plt.figure(figsize = (16, 7.5))
    plt.bar(x, y)
    plt.title("Runtime")
    plt.xlabel("Transformation")
    plt.xticks(rotation = 90)
    plt.ylabel("Runtime (s)")
    plt.subplots_adjust(left = 0.125, bottom = 0.28, right = 0.9, top = 0.88, wspace = 0.2, hspace = 0.2)
    plt.savefig(sample_id + "_runtime.png")

def plot_memory(transformations : list, sample_id : str):
    x = []
    y = []
    for i in range(1, len(transformations)):
        x.append(transformations[i].name)
        y.append(transformations[i].memory)
    plt.figure(figsize = (16, 7.5))
    plt.bar(x, y)
    plt.title("Peak Memory Usage")
    plt.xlabel("Transformation")
    plt.xticks(rotation = 90)
    plt.ylabel("Peak Memory Usage (MB)")
    plt.subplots_adjust(left = 0.125, bottom = 0.28, right = 0.9, top = 0.88, wspace = 0.2, hspace = 0.2)
    plt.savefig(sample_id + "_memory.png")

def plot_gene_mean_var(transformations : list, sample_id : str):
    fig = plt.figure(figsize = (16, 7.5))
    for t in range(0, len(transformations)):
        ax = fig.add_subplot(4, 5, t + 1)
        ax.scatter(transformations[t].matrix.mean(0), transformations[t].matrix.var(0), s = 2)
        ax.set_xlabel("Gene Mean")
        ax.set_ylabel("Gene Variance")
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        ax.legend([transformations[t].name], loc = 1, fontsize = "6")
    plt.subplots_adjust(left = 0.125, bottom = 0.03, right = 0.9, top = 0.98, wspace = 0.2, hspace = 0.2)
    plt.savefig(sample_id + "_gene_mean_var.png")

def plot_marker_genes(transformations : list, genes : dict, sample_id : str):
    for gene in genes.keys():
        fig = plt.figure(figsize = (16, 7.5))
        for t in range(0, len(transformations)):
            ax = fig.add_subplot(4, 5, t + 1)
            ax.hist(transformations[t].matrix.loc[:, gene], bins = 20)
            ax.set_xlabel("Gene Mean")
            ax.set_ylabel("No. Spots")
            ax.legend([transformations[t].name], loc = 1, fontsize = "6")
        plt.subplots_adjust(left = 0.05, bottom = 0.06, right = 0.99, top = 0.94, wspace = 0.34, hspace = 0.4)
        plt.suptitle(gene + "Distributions")
        plt.savefig(sample_id + "_" + gene + "_distributions.png")

def plot_hvg_similarity(transformations : list, name : list, sample_id : str):
    sims = pandas.DataFrame(0.0, index = name, columns = name)
    for i in transformations:
        i_hvg = i.hvgenes
        for j in transformations:
            j_hvg = j.hvgenes
            intersection = len(set(i_hvg).intersection(set(j_hvg)))
            union = len(set(i_hvg).union(set(j_hvg)))
            jaccard = intersection / union
            sims.at[i.name, j.name] = jaccard
    plt.figure(figsize = (16, 7.5))
    seaborn.heatmap(sims, cmap = "Blues", cbar = True, xticklabels = True, yticklabels = True)
    plt.subplots_adjust(left = 0.13, bottom = 0.25, right = 1.0, top = 0.98, wspace = 0.2, hspace = 0.2)
    plt.savefig(sample_id + "_hvg_similarity.png")

def plot_pathway_enrichment(enrichment : pandas.DataFrame, sample_id : str):
    plt.figure(figsize = (16, 7.5))
    seaborn.heatmap(enrichment, cmap = "Blues", cbar = True, xticklabels = True, yticklabels = True)
    plt.xlabel("Transformation")
    plt.ylabel("Marker Pathway")
    plt.subplots_adjust(left = 0.13, bottom = 0.25, right = 1.0, top = 0.98, wspace = 0.2, hspace = 0.2)
    plt.savefig(sample_id + "_pathway_enrichment.png")

def plot_sf_correspondence(t_one : list, t_two : list, sample_id : str):
    fig = plt.figure(figsize = (16, 7.5))
    for t in range(0, len(t_one)):
        ax = fig.add_subplot(4, 5, t + 1)
        x = (t_one[t].matrix.sum(1) / (t_one[t].matrix.sum().sum() / len(t_one[t].matrix.index))).to_list()
        x = scipy.stats.zscore(x)
        y = (t_two[t].matrix.sum(1) / (t_two[t].matrix.sum().sum() / len(t_two[t].matrix.index))).to_list()
        y = scipy.stats.zscore(y)
        ax.scatter(x, y, s = 5)
        ax.axline((0, 0), slope = 1, color = "red")
        ax.set_xlabel("Full Gene Panel")
        ax.set_ylabel("5000 Gene Panel")
        ax.legend([t_one[t].name], loc = 1, fontsize = "6")
    plt.subplots_adjust(left = 0.05, bottom = 0.06, right = 0.99, top = 0.99, wspace = 0.4, hspace = 0.4)
    plt.savefig(sample_id + "_sf_correspondence.png")

def plot_count_correspondence(t_one : list, t_two : list, genes : dict, sample_id : str):
    for gene in genes.keys():
        fig = plt.figure(figsize = (16, 7.5))
        for t in range(0, len(t_one)):
            ax = fig.add_subplot(4, 5, t + 1)
            x = t_one[t].matrix.loc[:, gene].to_list()
            y = t_two[t].matrix.loc[:, gene].to_list()
            ax.scatter(x, y, s = 5)
            ax.axline((0, 0), slope = 1, color = "red")
            ax.set_xlabel("Full Gene Panel")
            ax.set_ylabel("5000 Gene Panel")
            ax.legend([t_one[t].name], loc = 1, fontsize = "6")
        plt.suptitle(gene)
        plt.subplots_adjust(left = 0.05, bottom = 0.06, right = 0.99, top = 0.94, wspace = 0.5, hspace = 0.4)
        plt.savefig(sample_id + "_" + gene + "_count_correspondence.png")

def plot_corr_matrix(transformations : list, sample_id : str):
    fig = plt.figure(figsize = (16, 7.5))
    for t in range(0, len(transformations)):
        ax = fig.add_subplot(4, 5, t + 1)
        plt.title(transformations[t].name)
        seaborn.heatmap(transformations[t].matrix.corr(), ax = ax, cmap = "Blues", cbar = True, xticklabels = False, yticklabels = False)
    plt.subplots_adjust(left = 0.02, bottom = 0.03, right = 0.97, top = 0.94, wspace = 0.25, hspace = 0.3)
    plt.savefig(sample_id + "_corr_matrix.png")

def plot_total_rank(t_one : list, t_two : list, name : list, genes : dict, weight_matrix : pandas.DataFrame, enrichment : pandas.DataFrame, sample_id : str) -> pandas.DataFrame:
    runtime = pandas.Series(ranker.rank_runtime(t_one), index = name, name = "Runtime")
    memory = pandas.Series(ranker.rank_memory(t_one), index = name, name = "Memory")
    gmv = pandas.Series(ranker.rank_gene_mean_var(t_one), index = name, name = "Gene Mean-Variance")
    mgd = pandas.Series(ranker.rank_marker_genes(t_one, name, genes, weight_matrix), index = name, name = "Marker Gene Distribution")
    hvg = pandas.Series(ranker.rank_hvg_similarity(t_one), index = name, name = "HVG Similarity")
    enr = pandas.Series(ranker.rank_pathway_enrichment(enrichment), index = name, name = "Pathway Enrichment")
    sf = pandas.Series(ranker.rank_sf_correspondence(t_one, t_two), index = name, name = "Size Factor Correspondence")
    nc = pandas.Series(ranker.rank_count_correspondence(t_one, t_two, genes), index = name, name = "Count Correspondence")
    cm = pandas.Series(ranker.rank_corr_matrix(t_one), index = name, name = "Gene-Gene Correlation")
    rank = pandas.concat([runtime, memory, gmv, mgd, hvg, enr, sf, nc, cm], axis = 1)
    rank_sum = pandas.Series(rank.sum(1), index = name, name = "Sum")
    rank = pandas.concat([rank, rank_sum], axis = 1)
    rank.sort_values(by = "Sum", axis = 0, ascending = False, inplace = True)
    trans = rank.index.to_list()
    trans.remove("y")
    analyses = rank.columns.to_list()
    analyses.remove("Sum")
    plt.figure(figsize = (16, 7.5))
    seaborn.heatmap(rank.loc[trans, analyses], cmap = "Blues", cbar = True, xticklabels = True, yticklabels = True)
    plt.ylabel("Transformation")
    plt.subplots_adjust(left = 0.15, bottom = 0.28, right = 1.0, top = 0.98, wspace = 0.2, hspace = 0.2)
    plt.savefig(sample_id + "_total_rank.png")