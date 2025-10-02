import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys
from sklearn.metrics.cluster import adjusted_rand_score
import STAGATE_pyG
import tensorflow as tf
import json
from matplotlib.image import imread
from pathlib import Path

def chech_fileexist_stagate(args, filepath_output, slice_num):
    normalization = 'original_stagate'
    output_path = filepath_output + normalization +'/'
    output_path_cp = output_path + args.technology + '_' + args.dataset + '_' + slice_num + '_' 
    output_path_cp_cp = output_path_cp + 'run5_stagate_ari_score.txt'
    return os.path.exists(output_path_cp_cp)

def run_stagate_single(adata, output_path_cp, k_cutoff, alpha, n_clusters, seed):
    
    # stagate training
    STAGATE_pyG.Cal_Spatial_Net(adata, k_cutoff=k_cutoff, model='KNN')
    STAGATE_pyG.Stats_Spatial_Net(adata)
    tf.compat.v1.disable_eager_execution() # force eager
    adata = STAGATE_pyG.train_STAGATE(adata, random_seed=seed)
    
    # clustering
    adata = STAGATE_pyG.mclust_R(adata, used_obsm='STAGATE', num_cluster=n_clusters)
    obs_df = adata.obs.dropna()
    
    # calculate ari score
    obs_df = adata.obs.dropna()
    ARI = adjusted_rand_score(obs_df['mclust'], obs_df['ground_truth'])
    
    # save ari score
    try:
        ari_filename = '%sstagate_ari_score.txt' % (output_path_cp)
        with open(ari_filename, 'w') as f:
            f.write(f'ARI Score: {ARI}\n')
    except Exception as e:
        print(f"Failed to save ARI score without refinement: {e}")
    
    # save clustering result
    results = adata.obs['mclust']
    results.to_csv(output_path_cp + 'stagate_domain.csv', sep=',', header=True, index=True)
    
    # save embd
    pd.DataFrame(adata.obsm['STAGATE'], index=adata.obs_names).to_csv(output_path_cp + "stagate_latent_z.csv", index_label="barcode")

    # save umap
    sc.pp.neighbors(adata, use_rep='STAGATE')
    sc.tl.umap(adata)
    plt.rcParams["figure.figsize"] = (3, 3)
    sc.pl.umap(adata, color=["mclust", "ground_truth"], title=['STAGATE (ARI=%.2f)'%ARI, "Ground Truth"])
    try:
        plt.savefig(output_path_cp + "stagate_umap_predicted_labels.png", dpi=300)
    except Exception as e:
        print(f"Failed to save UMAP predicted labels plot: {e}")    

    # save physical plot
    sc.pl.spatial(adata,
                  img_key="hires",
                  color=["ground_truth", "mclust"],
                  title=["Ground truth", "ARI=%.4f"%ARI],
                  show=False)
    try:
        plt.savefig(output_path_cp + "stagate_spatial_clustering_result.png", dpi=300)
    except Exception as e:
        print(f"Failed to save spatial clustering result plot: {e}")

def run_stagate(args, file_dic, file_dic_gt, slice_num, filepath_output, hvg_type, n_top_genes = 3000, k_cutoff = 6, alpha=0):
    os.environ['R_HOME'] = '/oscar/rt/9.2/software/0.20-generic/0.20.1/opt/spack/linux-rhel9-x86_64_v3/gcc-11.3.1/r-4.3.1-lmofgb4ggmztfpsknzgazhhiwaua5ocd/rlib/R'
   
    
    filepath_norm = file_dic + args.dataset + '_' + slice_num + '/norm_counts/'
    norm_list = [f for f in os.listdir(filepath_norm) if os.path.isfile(os.path.join(filepath_norm, f))]
    norm_list = [os.path.splitext(f)[0] for f in norm_list if f.endswith('.h5ad')]
    norm_list = [n for n in norm_list if n not in ('scanpy_seurat', 'scanpy_weinreb')]
    norm_list.append("original_stagate")
    
    for normalization in norm_list:
        # define output path
        output_path = filepath_output + normalization +'/'

        # if not os.path.exists(check_path):
        if normalization == "original_stagate":
            filepath_input = file_dic + args.dataset + '_' + slice_num + '/norm_counts/raw.h5ad'
            adata = sc.read_h5ad(filepath_input)
        else:
            filepath_input = filepath_norm + normalization + '.h5ad'
            hvg_path = f"/oscar/data/yma16/Project/spTransform/1.Evaluations/10xVisium/DLPFC_{slice_num}/{hvg_type}_hvg_sdd/{normalization}.csv"
            adata = sc.read_h5ad(filepath_input)
            # subset hvg
            hvg_df = pd.read_csv(hvg_path,index_col=0)
            hvg_list = hvg_df.iloc[:, 0].tolist()
            adata = adata[:, adata.var_names.isin(hvg_list)].copy()
        spatial = pd.read_csv(file_dic + args.dataset + '_' + slice_num + '/coordinates.csv')
        adata.obsm["spatial"] = spatial.to_numpy()
        adata.var_names_make_unique()
        
        # read ground truth
        df_anno = pd.read_csv(file_dic_gt+'layers.csv')
        df_anno.set_index('Unnamed: 0', inplace=True)
        df_anno.index = df_anno.index.str.replace('\..*$', '', regex=True)
        df_anno_cp = df_anno[df_anno['sample_id'] == int(slice_num)]

        # add ground_truth
        if not adata.obs.index.equals(df_anno_cp.index):
            df_anno_cp = df_anno_cp.loc[adata.obs.index]
        adata.obs['ground_truth'] = df_anno_cp['layer_guess_reordered']
    
        # filter out NA nodes
        adata = adata[~pd.isnull(adata.obs['ground_truth'])]
        
        n_clusters = len(set(adata.obs['ground_truth'].unique()))
        
        #Preprocess
        if normalization == "original_stagate":
            sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n_top_genes)
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
        else:
            adata.var['highly_variable'] = True
            
        print(normalization, flush=True)
        print("hvg selection finished", flush=True)
        
        if not os.path.exists(output_path):
            os.makedirs(output_path)
            
        output_path_cp = output_path + args.technology + '_' + args.dataset + '_' + slice_num + '_' 
        
        if normalization == "original_stagate":
            highly_variable_genes = adata.var['highly_variable']
            highly_variable_genes_df = highly_variable_genes.to_frame()
            highly_variable_genes_df.to_csv(output_path_cp + 'highly_variable_genes.csv')
        
            print("hvg saved", flush=True)
        
        # add image info 
        path = Path(file_dic_gt+str(slice_num))
        
        library_id = str(slice_num)
        adata.uns["spatial"] = dict()
        adata.uns["spatial"][library_id] = dict()
        tissue_positions_file = (
            path / "spatial/tissue_positions.csv"
            if (path / "spatial/tissue_positions.csv").exists()
            else path / "spatial/tissue_positions_list.csv"
        )
        files = dict(
            tissue_positions_file=tissue_positions_file,
            scalefactors_json_file=path / "spatial/scalefactors_json.json",
            hires_image=path / "spatial/tissue_hires_image.png",
            lowres_image=path / "spatial/tissue_lowres_image.png",
        )
        
        adata.uns["spatial"][library_id]["images"] = dict()
        for res in ["hires", "lowres"]:
            try:
                adata.uns["spatial"][library_id]["images"][res] = imread(
                    str(files[f"{res}_image"])
                )
            except Exception:
                raise OSError(f"Could not find '{res}_image'")
        
        # read json scalefactors
        adata.uns["spatial"][library_id]["scalefactors"] = json.loads(
            files["scalefactors_json_file"].read_bytes()
        )
        
        # read coordinates
        positions = pd.read_csv(
            files["tissue_positions_file"],
            header=0 if tissue_positions_file.name == "tissue_positions.csv" else None,
            index_col=0,
        )
        positions.columns = [
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_col_in_fullres",
            "pxl_row_in_fullres",
        ]
        
        adata.obs = adata.obs.join(positions, how="left")
        
        adata.obsm["spatial"] = adata.obs[
            ["pxl_row_in_fullres", "pxl_col_in_fullres"]
        ].to_numpy()
        adata.obs.drop(
            columns=["pxl_row_in_fullres", "pxl_col_in_fullres"],
            inplace=True,
        )  
        
        print('images saved', flush=True)
        ##### Run
        for i in range(1, 6):
            output_path_cp2 = output_path_cp + 'run' + str(i) + '_'
            print('run' + str(i) +  'start', flush=True)
            try:
                run_stagate_single(adata, output_path_cp2, k_cutoff, alpha, n_clusters, i)
            except Exception as e:
                print('Error encountered at ' + normalization)
                print(f"An error occurred: {e}")
            print('run' + str(i) +  'finish', flush=True)