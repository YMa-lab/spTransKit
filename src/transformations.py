import pandas
import numpy
import reactome2py.analysis
import scanpy
import random
import argparse
import os
import shutil

import plotter

def handle_coordinates(file : str, technology : str) -> pandas.DataFrame:
    if technology == "10xVisium":
        coordinates = pandas.read_csv(file, index_col = 0, header = 0)
        coordinates = coordinates.loc[coordinates["in_tissue"] == 1, ["array_col", "array_row"]]
        coordinates.rename(columns = {"array_col" : "x", "array_row" : "y"}, inplace = True)
    if technology == "Slide-seq":
        coordinates = pandas.read_csv(file, index_col = 0, header = 0)
        coordinates = coordinates.loc[:, ["xcoord", "ycoord"]]
        coordinates.rename(columns = {"xcoord" : "x", "ycoord" : "y"}, inplace = True)
    if technology == "Merfish":
        coordinates = pandas.read_csv(file, index_col = 0, header = 0)
        coordinates = coordinates.loc[:, ["x", "y"]]
    if technology == "Stereo-seq":
        coordinates = scanpy.read_h5ad(file).to_df()
    return coordinates

def handle_matrix(file : str, coordinates : pandas.DataFrame, technology : str) -> pandas.DataFrame:
    if technology == "10xVisium":
        matrix = scanpy.read_10x_h5(file).to_df()
        matrix = matrix.loc[coordinates.index, :]   
    if technology == "Slide-seq":
        matrix = pandas.read_csv(file, index_col = 0, header = 0).T
    if technology == "Merfish":
        matrix = pandas.read_csv(file, index_col = 0, header = 0)
    if technology == "Stereo-seq":
        matrix = scanpy.read_h5ad(file).to_df()
   
    matrix = matrix.loc[:, (matrix != 0).any(axis = 0)]
    matrix = matrix.loc[:, ~matrix.columns.duplicated()]
    matrix = matrix.loc[(matrix != 0).any(axis = 1), :]
    coordinates = coordinates.loc[matrix.index, :]

    return matrix

def handle_marker_genes(file : str) -> dict:
    genes = pandas.read_csv(file, header = 0)
    return pandas.Series(data = genes.distribution.values, index = genes.gene).to_dict()

def handle_marker_pathways(file : str):
    return pandas.read_csv(file, header = None).loc[:, 0].to_list()

def transform_test(counts : pandas.DataFrame, num_hvg : int) -> list:
    from functions import raw, rawsize, cpm, normalisr
    one = raw.Raw(counts, num_hvg)
    two = rawsize.RawSize(counts, num_hvg)
    three = cpm.CPM(counts, num_hvg)
    four = normalisr.Normalisr(counts, num_hvg)
    return [one, two, three, four]

def transform(counts : pandas.DataFrame, num_hvg : int) -> list:
    from functions import raw, rawsize, cpm, shiftedlog, acosh, pseudoshiftedlog, cpmshiftedlog, shiftedlogsize, analyticpearson, scanpy_zheng, scanpy_weinreb, scanpy_seurat, scanpy_pearson_residual, deseq2, tmm, normalisr, psinorm
    one = raw.Raw(counts, num_hvg)
    two = rawsize.RawSize(counts, num_hvg)
    three = cpm.CPM(counts, num_hvg)
    four = shiftedlog.ShiftedLog(counts, num_hvg)
    five = acosh.Acosh(counts, num_hvg)
    six = pseudoshiftedlog.PseudoShiftedLog(counts, num_hvg)
    seven = cpmshiftedlog.CPMShiftedLog(counts, num_hvg)
    eight = shiftedlogsize.ShiftedLogSize(counts, num_hvg)
    nine = analyticpearson.PearsonResidual(counts, num_hvg)
    ten = scanpy_zheng.ScanpyZheng(counts, num_hvg)
    eleven = scanpy_weinreb.ScanpyWeinreb(counts, num_hvg)
    twelve = scanpy_seurat.ScanpySeurat(counts, num_hvg)
    thirteen = scanpy_pearson_residual.ScanpyPearsonResidual(counts, num_hvg)
    fourteen = deseq2.DESeq2(counts, num_hvg)
    fifteen = tmm.TMM(counts, num_hvg)
    sixteen = normalisr.Normalisr(counts, num_hvg)
    seventeen = psinorm.PsiNorm(counts, num_hvg)
    return [one, two, three, four, five, six, seven, eight, nine, ten, eleven, twelve, thirteen, fourteen, fifteen, sixteen, seventeen]

def names(transformations : list) -> list:
    name = []
    for x in transformations:
        name.append(x.name)
    return name

def gene_enrichment(transformations : list, name : list, pwys : list) -> pandas.DataFrame:
    y = pandas.DataFrame(1.0, index = pwys, columns = name)
    for t in transformations:
        hvg = ""
        for gene in t.hvgenes:
            hvg = hvg + gene + ","
        result = reactome2py.analysis.identifiers(ids = hvg, p_value = "0.05")
        token = result["summary"]["token"]
        marker_pathways = ""
        for pwy in pwys:
            marker_pathways = marker_pathways + pwy + ","
        pathways = reactome2py.analysis.token_pathways_result(token, marker_pathways)
        pathways_series = pandas.Series(1.0, index = pwys)
        for pwy in pathways:
            pathways_series[pwy["stId"]] = float(pwy["entities"]["pValue"])
        y[t.name] = pathways_series
    return -pandas.DataFrame(numpy.log10(y))

def transform_again_test(counts: pandas.DataFrame, num_hvg : int, subpanel : list) -> list:
    from functions import raw, rawsize, cpm, normalisr
    sub_counts = counts.loc[:, subpanel]
    sub_counts = sub_counts.loc[(sub_counts != 0).any(axis = 1), :]
    sub_counts = sub_counts.loc[(sub_counts != 0.0).any(axis = 1), :]
    one = raw.Raw(sub_counts, num_hvg)
    two = rawsize.RawSize(sub_counts, num_hvg)
    three = cpm.CPM(sub_counts, num_hvg)
    four = normalisr.Normalisr(sub_counts, num_hvg)
    return [one, two, three, four]

def transform_again(counts : pandas.DataFrame, num_hvg : int, genes : list) -> list:
    from functions import raw, rawsize, cpm, shiftedlog, acosh, pseudoshiftedlog, cpmshiftedlog, shiftedlogsize, analyticpearson, scanpy_zheng, scanpy_weinreb, scanpy_seurat, scanpy_pearson_residual, deseq2, tmm, normalisr, psinorm
    sub_counts = counts.loc[:, genes]
    sub_counts = sub_counts.loc[(sub_counts != 0).any(axis = 1), :]
    sub_counts = sub_counts.loc[(sub_counts != 0.0).any(axis = 1), :]
    one = raw.Raw(sub_counts, num_hvg)
    two = rawsize.RawSize(sub_counts, num_hvg)
    three = cpm.CPM(sub_counts, num_hvg)
    four = shiftedlog.ShiftedLog(sub_counts, num_hvg)
    five = acosh.Acosh(sub_counts, num_hvg)
    six = pseudoshiftedlog.PseudoShiftedLog(sub_counts, num_hvg)
    seven = cpmshiftedlog.CPMShiftedLog(sub_counts, num_hvg)
    eight = shiftedlogsize.ShiftedLogSize(sub_counts, num_hvg)
    nine = analyticpearson.PearsonResidual(sub_counts, num_hvg)
    ten = scanpy_zheng.ScanpyZheng(sub_counts, num_hvg)
    eleven = scanpy_weinreb.ScanpyWeinreb(sub_counts, num_hvg)
    twelve = scanpy_seurat.ScanpySeurat(sub_counts, num_hvg)
    thirteen = scanpy_pearson_residual.ScanpyPearsonResidual(sub_counts, num_hvg)
    fourteen = deseq2.DESeq2(sub_counts, num_hvg)
    fifteen = tmm.TMM(sub_counts, num_hvg)
    sixteen = normalisr.Normalisr(sub_counts, num_hvg)
    seventeen = psinorm.PsiNorm(sub_counts, num_hvg)
    return [one, two, three, four, five, six, seven, eight, nine, ten, eleven, twelve, thirteen, fourteen, fifteen, sixteen, seventeen]

def weight_matrix(coordinates : pandas.DataFrame) -> pandas.DataFrame:
    x = pandas.DataFrame(0, index = coordinates.index, columns = coordinates.index)
    x = x.add(coordinates.loc[:, "x"], axis = 0).sub(coordinates.loc[:, "x"], axis = 1) ** 2
    y = pandas.DataFrame(0, index = coordinates.index, columns = coordinates.index)
    y = y.add(coordinates.loc[:, "y"], axis = 0).sub(coordinates.loc[:, "y"], axis = 1) ** 2
    weight_matrix = pandas.DataFrame(numpy.exp(-(x.add(y) ** 0.5)))
    return weight_matrix

def create_subpanel(matrix : pandas.DataFrame, genes : dict) -> list:
    subpanel = random.sample(matrix.columns.to_list(), int(len(matrix.columns.to_list()) / 4))
    for gene in genes.keys():
        if gene not in subpanel:
            subpanel.pop(0)
            subpanel.append(gene)
    return subpanel

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_id", type = str, help = "sample identifier")
    parser.add_argument("--technology", type = str, help = "type of ST technology")
    parser.add_argument("--counts", type = str, help = "file with gene counts")
    parser.add_argument("--coordinates", type = str, help = "file with spot coordinates")
    parser.add_argument("--genes", type = str, help = "file with marker genes")
    parser.add_argument("--pathways", type = str, help = "file with marker pathways")
    parser.add_argument("--num_hvg", type = int, help = "number of highly variable genes to be selected")
    params = parser.parse_args()

    coordinates = handle_coordinates(params.coordinates, params.technology)
    matrix = handle_matrix(params.counts, coordinates, params.technology)
    genes = handle_marker_genes(params.genes)
    pathways = handle_marker_pathways(params.pathways)

    t_one = transform(matrix, params.num_hvg)
    name = names(t_one)

    subpanel = create_subpanel(matrix, genes)
    t_two = transform_again(matrix, params.num_hvg, subpanel)

    enrichment = gene_enrichment(t_one, name, pathways)

    wm = weight_matrix(coordinates)

    dir = params.technology + "/" + params.sample_id

    if os.path.exists(dir):
        shutil.rmtree(dir)
    os.makedirs(dir)

    plotter.plot_runtime(t_one, dir)
    plotter.plot_memory(t_one, dir)
    plotter.plot_gene_mean_var(t_one, dir)
    plotter.plot_marker_genes(t_one, genes, dir)
    plotter.plot_hvg_similarity(t_one, name, dir)
    plotter.plot_pathway_enrichment(enrichment, dir)
    plotter.plot_sf_correspondence(t_one, t_two, dir)
    plotter.plot_count_correspondence(t_one, t_two, genes, dir)
    plotter.plot_corr_matrix(t_one, dir)
    plotter.plot_total_rank(t_one, t_two, name, genes, wm, enrichment, dir)
