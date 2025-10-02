import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import SpaGCN as spg
import cv2
from sklearn import metrics
from pathlib import Path
from matplotlib.image import imread
import json

def chech_fileexist_spagcn(args, filepath_output, slice_num):
    normalization = 'original_spagcn'
    output_path = filepath_output + normalization +'/'
    output_path_cp = output_path + args.technology + '_' + args.dataset + '_' + slice_num + '_' 
    output_path_cp_cp = output_path_cp + 'run5_spagcn_ari_score_wrefine.txt'
    return os.path.exists(output_path_cp_cp)

def run_spagcn_single(adata, adj, adj_2d, l, n_clusters, seed, output_path_cp):
    r_seed=t_seed=n_seed=seed
    res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)
    
    clf=spg.SpaGCN()
    clf.set_l(l)
    #Set seed
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    #Run
    clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
    embed, y_pred, prob=clf.predict()

    print(f"embed shape{embed.shape}")

    adata.obs["pred"]= y_pred
    adata.obs["pred"]=adata.obs["pred"].astype('category')
    
    refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
    adata.obs["refined_pred"]=refined_pred
    adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
    
    adata.var.columns = adata.var.columns.astype(str)
    
    # calculate metric ARI
    ARI_worefine = metrics.adjusted_rand_score(adata.obs['pred'], adata.obs['ground_truth'])
    adata.uns['ARI_worefine'] = ARI_worefine
    
    # plotting spatial clustering result
    sc.pl.spatial(adata,
                  img_key="hires",
                  color=["ground_truth", "pred"],
                  title=["Ground truth", "ARI without refinement=%.4f"%ARI_worefine],
                  show=False)
    
    try:
        plt.savefig(output_path_cp + "spagcn_spatial_clustering_result_worefine.png", dpi=300)
        plt.close()
    except Exception as e:
        print(f"Failed to save spatial clustering result worefine plot: {e}")

    # calculate metric ARI
    ARI_wrefine = metrics.adjusted_rand_score(adata.obs['refined_pred'], adata.obs['ground_truth'])
    adata.uns['ARI_wrefine'] = ARI_wrefine
    
    # plotting spatial clustering result
    sc.pl.spatial(adata,
                  img_key="hires",
                  color=["ground_truth", "refined_pred"],
                  title=["Ground truth", "ARI with refinement=%.4f"%ARI_wrefine],
                  show=False)
    
    # save embed
    if isinstance(embed, torch.Tensor):
        embed = embed.detach().cpu().numpy()
    else:
        embed = np.asarray(embed)
    pd.DataFrame(embed, index=adata.obs_names).to_csv(output_path_cp + "spagcn_latent_z.csv", index_label="barcode")

    try:
        plt.savefig(output_path_cp + "spagcn_spatial_clustering_result_wrefine.png", dpi=300)
        plt.close()
    except Exception as e:
        print(f"Failed to save spatial clustering result wrefine plot: {e}")

    try:
        results_worefine = adata.obs['pred']
        results_worefine.to_csv(output_path_cp + 'spagcn_domain_worefine.csv', sep=',', header=True, index=True)
    except Exception as e:
        print(f"Failed to save domain without refinement results: {e}")
        
    try:
        results_wrefine = adata.obs['refined_pred']
        results_wrefine.to_csv(output_path_cp + 'spagcn_domain_wrefine.csv', sep=',', header=True, index=True)
    except Exception as e:
        print(f"Failed to save domain without refinement results: {e}")

    try:
        ari_filename = '%sspagcn_ari_score_worefine.txt' % (output_path_cp)
        with open(ari_filename, 'w') as f:
            f.write(f'ARI Score without refinement: {ARI_worefine}\n')
    except Exception as e:
        print(f"Failed to save ARI score without refinement: {e}")
        
    try:
        ari_filename = '%sspagcn_ari_score_wrefine.txt' % (output_path_cp)
        with open(ari_filename, 'w') as f:
            f.write(f'ARI Score with refinement: {ARI_wrefine}\n')
    except Exception as e:
        print(f"Failed to save ARI score with refinement: {e}")

    return None

def run_spagcn(args, file_dic, file_dic_gt, slice_num, filepath_output, hvg_type, s=1, b=49, p=0.5):

    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    
    os.environ['R_HOME'] = '/oscar/rt/9.2/software/0.20-generic/0.20.1/opt/spack/linux-rhel9-x86_64_v3/gcc-11.3.1/r-4.3.1-lmofgb4ggmztfpsknzgazhhiwaua5ocd/rlib/R' 
    
    # read data
    filepath_norm = file_dic + args.dataset + '_' + slice_num + '/norm_counts/'
    norm_list = [f for f in os.listdir(filepath_norm) if os.path.isfile(os.path.join(filepath_norm, f))]
    norm_list = [os.path.splitext(f)[0] for f in norm_list if f.endswith('.h5ad')]
    norm_list = [n for n in norm_list if n not in ('scanpy_seurat', 'scanpy_weinreb')]
    norm_list.append("original_spagcn")
    
    for normalization in norm_list:
        print(normalization, flush=True)
        if normalization == "original_spagcn":
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
        adata.var_names_make_unique()
        
        # load pixel_pos and array_pos
        pixel_pos = pd.read_csv(file_dic_gt+'pixel_pos.csv', index_col=0)
        array_pos = pd.read_csv(file_dic_gt+'array_pos.csv', index_col=0)
        array_pos.index = array_pos.index.str.split('.').str[0]
        pos_all = pd.concat([pixel_pos, array_pos], axis=1)
        pos_slice = pos_all[pos_all['sample_id']==int(slice_num)]
        
        adata.obs["x_array"]=pos_slice['array_row']
        adata.obs["y_array"]=pos_slice['array_col']
        adata.obs["x_pixel"]=pos_slice['pxl_row_in_fullres']
        adata.obs["y_pixel"]=pos_slice['pxl_col_in_fullres']
            
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
    
        # load image
        img=cv2.imread(f"{file_dic_gt}{slice_num}/{slice_num}_full_image.tif")
        
        #Set coordinates
        x_array=adata.obs["x_array"].tolist()
        y_array=adata.obs["y_array"].tolist()
        x_pixel=adata.obs["x_pixel"].tolist()
        y_pixel=adata.obs["y_pixel"].tolist()
        
        # 
        img_path = f"{file_dic_gt}{slice_num}/{slice_num}_map.jpg"
        if not os.path.exists(img_path):
            
            #Test coordinates on the image
            img_new=img.copy()
            for i in range(len(x_pixel)):
                x=x_pixel[i]
                y=y_pixel[i]
                img_new[int(x-20):int(x+20), int(y-20):int(y+20),:]=0
            cv2.imwrite(img_path, img_new)
        
        print('all data loaded saved', flush=True)
        
        adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, image=img, beta=b, alpha=s, histology=True)
        adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
        
        spg.prefilter_genes(adata,min_cells=3) 
        spg.prefilter_specialgenes(adata)
        print("num of cells after filtering = " + str(adata.shape[0]))
        print("num of genes after filtering = " + str(adata.shape[1]))
        if normalization == "original_spagcn":
            #Normalize and take log for UMI
            sc.pp.normalize_per_cell(adata)
            sc.pp.log1p(adata)   
           
        l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)
        
        for i in range(1, 6):
            output_path_cp2 = output_path_cp + 'run' + str(i) + '_'
            print('run' + str(i) +  'start', flush=True)
            try:
                run_spagcn_single(adata, adj, adj_2d, l, n_clusters, i, output_path_cp2)
            except Exception as e:
                print('Error encountered at ' + normalization)
                print(f"An error occurred: {e}")
            print('run' + str(i) +  'finish', flush=True)
        

