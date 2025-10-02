import os
import torch
import pandas as pd
import matplotlib
matplotlib.use('template')  # or another valid backend from the list
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp
from GraphST import GraphST
from GraphST.utils import clustering
from pathlib import Path
from matplotlib.image import imread
import json
import matplotlib.pyplot as plt

def chech_fileexist(args, filepath_output, slice_num):
    normalization = 'original_graphst'
    output_path = filepath_output + normalization +'/'
    output_path_cp = output_path + args.technology + '_' + args.dataset + '_' + slice_num + '_' 
    output_path_cp_cp = output_path_cp + 'run5_graphst_ari_score.txt'
    return os.path.exists(output_path_cp_cp)
    
def run_graphst_single(adata, output_path_cp, seed, types, device, n_clusters, tool, radius):
    # define model
    if types == 'original':
        model = GraphST.GraphST(adata, device=device, random_seed = seed, types='original')
    else:
        model = GraphST.GraphST(adata, device=device, random_seed = seed)
    
    
    # train model
    adata = model.train()
    
    # clustering
    if tool == 'mclust':
        clustering(adata, n_clusters, radius=radius, method=tool, refinement=True)
    elif tool in ['leiden', 'louvain']:
        clustering(adata, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False)
    
    
    # calculate metric ARI
    ARI = metrics.adjusted_rand_score(adata.obs['domain'], adata.obs['ground_truth'])
    adata.uns['ARI'] = ARI
    
    # plotting spatial clustering result
    sc.pl.spatial(adata,
                  img_key="hires",
                  color=["ground_truth", "domain"],
                  title=["Ground truth", "ARI=%.4f"%ARI],
                  show=False)
    
    # save embedding
    pd.DataFrame(adata.obsm['emb'], index=adata.obs_names).to_csv(f"{output_path_cp}graphst_embedding.csv", index_label="barcode")

    # plt.savefig(output_path_cp + "graphst_spatial_clustering_result.png", dpi=300)
    try:
        plt.savefig(output_path_cp + "graphst_spatial_clustering_result.png", dpi=300)
    except Exception as e:
        print(f"Failed to save spatial clustering result plot: {e}")

    # plotting predicted labels by UMAP
    sc.pp.neighbors(adata, use_rep='emb_pca', n_neighbors=10)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color='domain', title=['Predicted labels'], show=False)
    # plt.savefig(output_path_cp + "graphst_umap_predicted_labels.png", dpi=300)

    try:
        plt.savefig(output_path_cp + "graphst_umap_predicted_labels.png", dpi=300)
    except Exception as e:
        print(f"Failed to save UMAP predicted labels plot: {e}")
        

    try:
        results = adata.obs['domain']
        results.to_csv(output_path_cp + 'graphst_domain.csv', sep=',', header=True, index=True)
    except Exception as e:
        print(f"Failed to save domain results: {e}")


    try:
        ari_filename = '%sgraphst_ari_score.txt' % (output_path_cp)
        with open(ari_filename, 'w') as f:
            f.write(f'ARI Score: {ARI}\n')
    except Exception as e:
        print(f"Failed to save ARI score: {e}")

    return None

def run_graphst(args, file_dic, file_dic_gt, slice_num, filepath_output, hvg_type, radius = 50, tool = 'mclust'):

    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    
    os.environ['R_HOME'] = '/oscar/rt/9.2/software/0.20-generic/0.20.1/opt/spack/linux-rhel9-x86_64_v3/gcc-11.3.1/r-4.3.1-lmofgb4ggmztfpsknzgazhhiwaua5ocd/rlib/R' 
    
    # read data
    filepath_norm = file_dic + args.dataset + '_' + slice_num + '/norm_counts/'
    norm_list = [f for f in os.listdir(filepath_norm) if os.path.isfile(os.path.join(filepath_norm, f))]
    norm_list = [os.path.splitext(f)[0] for f in norm_list if f.endswith('.h5ad')]
    norm_list.append("original_graphst")
    norm_list = [n for n in norm_list if n not in ('scanpy_seurat', 'scanpy_weinreb')]
    print(norm_list)
    
    for normalization in norm_list:
        print(normalization, flush=True)
        # define output path
        output_path = filepath_output + normalization +'/'
        fname = f"10xVisium_DLPFC_{slice_num}_run5_graphst_ari_score.txt"
        check_path = os.path.join(output_path, fname)

        if not os.path.exists(check_path):
            if normalization == "original_graphst":
                filepath_input = file_dic + args.dataset + '_' + slice_num + '/norm_counts/raw.h5ad'
                adata = sc.read_h5ad(filepath_input)
            else:
                filepath_input = filepath_norm + normalization + '.h5ad'
                hvg_path = f"/oscar/data/yma16/Project/spTransform/1.Evaluations/10xVisium/DLPFC_{slice_num}/{hvg_type}_hvg_sdd/{normalization}.csv"
                adata = sc.read_h5ad(filepath_input)
                # subset hvg
                hvg_df = pd.read_csv(hvg_path,index_col=0)
                hvg_list = hvg_df.iloc[:, 0].tolist()
                print(f'hvg_len is {len(hvg_list)}')
                adata = adata[:, adata.var_names.isin(hvg_list)].copy()

            spatial = pd.read_csv(file_dic + args.dataset + '_' + slice_num + '/coordinates.csv')
            adata.obsm["spatial"] = spatial.to_numpy()
            print('count shape is '+ str(adata.shape))
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
            
            # define output path
            output_path = filepath_output + normalization +'/'
            
            if not os.path.exists(output_path):
                os.makedirs(output_path)
                
            output_path_cp = output_path + args.technology + '_' + args.dataset + '_' + slice_num + '_'
            
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
            
            print('all data loaded saved', flush=True)
            
            for i in range(1, 7):
                output_path_cp2 = output_path_cp + 'run' + str(i) + '_'
                print('run' + str(i) +  'start', flush=True)
                try:
                    if normalization == "original_graphst":
                        run_graphst_single(adata, output_path_cp2, i, 'original', device, n_clusters, tool, radius)
                    else:
                        run_graphst_single(adata, output_path_cp2, i, 'new', device, n_clusters, tool, radius)
                except Exception as e:
                    print('Error encountered at ' + normalization)
                    print(f"An error occurred: {e}")
                print('run' + str(i) +  'finish', flush=True)
        else:
            print(f"Skipping slice {slice_num} normalization {normalization}, output already exists.")
