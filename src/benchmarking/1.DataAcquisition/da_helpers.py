import pandas as pd
import scanpy as sc

def handle_coordinates(file : str, technology : str) -> pd.DataFrame:
    # Test dataset
    if technology == "test":
        coordinates = pd.read_csv(file, index_col = 0, header = 0)
    
    # 10x Visium datasets
    if technology == "10xVisium":
        coordinates = pd.read_csv(file, index_col = 0, header = 0)
        if "array_col" in coordinates.columns:
            coordinates = coordinates.loc[coordinates["in_tissue"] == 1, ["array_col", "array_row"]]
        else:
            coordinates = pd.read_csv(file, index_col = 0, header = None)
            coordinates.columns = ["in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
            coordinates = coordinates.loc[coordinates["in_tissue"] == 1, ["array_col", "array_row"]]
        coordinates.rename(columns = {"array_col" : "x", "array_row" : "y"}, inplace = True)
    
    # Open-ST datasets
    if technology == "Open-ST":
        coordinates = pd.read_csv(file, index_col = 0, header = 0)
        coordinates.index = coordinates.index.astype(str)

    # Slide-seq and Slide-seqV2 datasets
    if technology == "Slide-seq" or technology == "Slide-seqV2":
        coordinates = pd.read_csv(file, index_col = 0, header = 0)
        if 1 in coordinates.index.to_list():
            coordinates.set_index("barcode", inplace = True)
        if "x" in coordinates.columns.to_list() and "y" in coordinates.columns.to_list():
            coordinates = coordinates.loc[:, ["x", "y"]]
        else:
            coordinates = coordinates.loc[:, ["xcoord", "ycoord"]]
            coordinates.rename(columns = {"xcoord" : "x", "ycoord" : "y"}, inplace = True)

    # Slide-Tag datasets
    if technology == "Slide-Tag":
        coordinates = pd.read_csv(file, index_col = 0, header = 0).drop(index = "TYPE")
        coordinates = coordinates.loc[:, ["X", "Y"]].astype("Float64")
        coordinates.rename(columns = {"X" : "x", "Y" : "y"}, inplace = True)
    
    # Stereo-seq datasets
    if technology == "Stereo-seq":
        coordinates = pd.read_csv(file, index_col = 0, header = 0)
        coordinates.index = coordinates.index.astype(str)

    # 10x Xenium and 10x Xenium Prime datasets
    if technology == "10x_Xenium" or technology == "10x_Xenium_Prime":
        coordinates = pd.read_csv(file, index_col = 0, header = 0)
        coordinates = coordinates.loc[:, ["x_centroid", "y_centroid"]]
        coordinates.rename(columns = {"x_centroid" : "x", "y_centroid" : "y"}, inplace = True)
        coordinates.index = coordinates.index.astype(str)

    # MERFISH datasets
    if technology == "Merfish":
        coordinates = pd.read_csv(file, index_col = 0, header = 0)
        coordinates = coordinates.loc[:, ["x", "y"]]
    
    # NanoString CosMx datasets
    if technology == "NanoString_CosMx":
        coordinates = pd.read_csv(file, index_col = None, header = 0)
        coordinates.index = coordinates["fov"].astype(str) + "_" + coordinates["cell_ID"].astype(str)
        coordinates = coordinates.loc[:, ["CenterX_global_px", "CenterY_global_px"]]
        coordinates.rename(columns = {"CenterX_global_px" : "x", "CenterY_global_px" : "y"}, inplace = True)

    # seqFISH datasets
    if technology == "seqFISH":
        coordinates = pd.read_csv(file, index_col = 0, header = 0)
    
    # seqFISH+ datasets
    if technology == "seqFISH+":
        coordinates = pd.read_csv(file, header = 0, index_col = 0)
        coordinates = coordinates.loc[:, ["sdimx", "sdimy"]]
        coordinates.rename(columns = {"sdimx" : "x", "sdimy" : "y"}, inplace = True)

    return coordinates

def handle_counts(file : str, coordinates : pd.DataFrame, technology : str) -> pd.DataFrame:
    # Test dataset
    if technology == "test":
        counts = pd.read_csv(file, index_col = 0, header = 0)
    
    # 10x Visium datasets
    if technology == "10xVisium":
        if ".h5" in file:
            counts = sc.read_10x_h5(file).to_df()
        elif ".csv" in file:
            counts = pd.read_csv(file, index_col = 0, header = 0)
        counts = counts.loc[coordinates.index, :]
        genes = counts.columns.to_list()
        counts.columns = [gene.replace(".", "-") for gene in genes]
    
    # Open-ST datasets
    if technology == "Open-ST":
        counts = sc.read_h5ad(file).to_df().astype(int)
        counts.index = counts.index.astype(str)

    # Slide-seq datasets
    if technology == "Slide-seq":
        counts = pd.read_csv(file, index_col = 0, header = 0)
        if 1 in counts.index.to_list():
            counts.set_index("barcode", inplace = True)
    
    # Slide-seqV2 datasets
    if technology == "Slide-seqV2":
        counts = pd.read_csv(file, index_col = 0, header = 0).T

    # Slide-Tag datasets
    if technology == "Slide-Tag":
        counts = pd.read_csv(file, index_col = 0, header = 0)
        cells = set(counts.index.to_list()).intersection(set(coordinates.index.to_list()))
        counts = counts.loc[list(cells), :]
        genes = counts.columns.to_list()
        counts.columns = [gene.replace(".", "-") for gene in genes]

    # Stereo-seq datasets
    if technology == "Stereo-seq":
        counts = sc.read_h5ad(file).to_df()
    
    # 10x Xenium and 10x Xenium Prime datasets
    if technology == "10x_Xenium" or technology == "10x_Xenium_Prime":
        counts = sc.read_10x_h5(file).to_df()
        counts.index = counts.index.astype(str)

    # MERFISH datasets
    if technology == "Merfish":
        counts = pd.read_csv(file, index_col = 0, header = 0)
    
    # NanoString CosMx datasets
    if technology == "NanoString_CosMx":
        counts = pd.read_csv(file, index_col = None, header = 0)
        counts.index = counts["fov"].astype(str) + "_" + counts["cell_ID"].astype(str)
        counts.drop(columns = ["fov", "cell_ID"], inplace = True)
        counts = counts.loc[coordinates.index.to_list(), :]

    # seqFISH datasets
    if technology == "seqFISH":
        counts = pd.read_csv(file, index_col = 0, header = 0)
    
    # seqFISH+ datasets
    if technology == "seqFISH+":
        counts = pd.read_csv(file, index_col = 0, header = 0)

    counts = counts.loc[:, (counts != 0).any(axis = 0)]
    counts = counts.loc[:, ~counts.columns.duplicated()]
    counts = counts.loc[(counts != 0).any(axis = 1), :]

    counts.index.name = "location"

    return counts

def check_dimensions(counts : pd.DataFrame, coordinates : pd.DataFrame):
    if counts.shape[0] != coordinates.shape[0]:
        raise Exception("Number of cells do not match between counts and coordinates!")
