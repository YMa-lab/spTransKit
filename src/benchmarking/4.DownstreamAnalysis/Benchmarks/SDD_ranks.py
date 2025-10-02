import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import os

def extract_ari_scores_all(technology, dataset, method, slice_num_list, suffix, hvg):
    ari_scores = {}
    
    # extract the save ARI scores by adjusting specific file format and suffix
    for slice_num in slice_num_list:
        base_path = f'/oscar/data/yma16/Project/spTransform/2.Downstream/results/{technology}/{dataset}_hvg_{hvg}_sdd/{slice_num}/{method}'

        for normalization in os.listdir(base_path):
            norm_dir = os.path.join(base_path, normalization)
            if os.path.isdir(norm_dir):
                if normalization not in ari_scores:
                    ari_scores[normalization] = {}
                if slice_num not in ari_scores[normalization]:
                    ari_scores[normalization][slice_num] = {}
                
                for i in range(1, 6):
                    if method == 'spagcn':
                        file_path = os.path.join(norm_dir, f"{technology}_{dataset}_{slice_num}_run{i}_{method}_ari_score_wrefine.txt")
                    elif method == 'bayesspace':
                        file_path = os.path.join(norm_dir, f"{technology}_{dataset}_run{i}_{method}_ari_score.txt")
                    else:
                        file_path = os.path.join(norm_dir, f"{technology}_{dataset}_{slice_num}_run{i}_{method}_ari_score.txt")
                    
                    if os.path.isfile(file_path):  
                        try:
                            with open(file_path, 'r') as file:
                                for line in file:
                                    text = line.strip()
                                    score_str = None
                                    if method == 'spagcn' and line.startswith("ARI Score with refinement:"):
                                        score_str = text.split(":", 1)[1]
                                    elif text.startswith("ARI Score:"):
                                        score_str = text.split(":", 1)[1]
                                    else:
                                        try:
                                            ari_score = float(text)
                                        except ValueError:
                                            # not a pure-number line—skip to next line
                                            continue
                                        else:
                                            # parsed successfully
                                            ari_scores[normalization][slice_num][f"run{i}"] = ari_score
                                            break
                                    if score_str is not None:
                                        try:
                                            ari_score = float(score_str.strip())
                                        except ValueError:
                                            # prefix found but not a valid float—skip
                                            continue
                                        ari_scores[normalization][slice_num][f"run{i}"] = ari_score
                                        break
                        except Exception as e:
                            print(f"Error reading file {file_path}: {e}")
                            ari_scores[normalization][slice_num][f"run{i}"] = 0.0
    
    return ari_scores

def extract_combined(hvg):
    technology='10xVisium'
    dataset='DLPFC'
    slice_num_list = ['151508', '151509', '151507', '151671', '151676', '151670', '151669', '151510', '151675', '151672', '151673', '151674']
    methods = ['stagate', 'spagcn', 'graphst','bayesspace']
    method_suffix = {'graphst': '_hvg', 'stagate': '_hvg', 'spagcn': '', 'bayesspace':'_hvg'}  # Define suffixes for each method
    method_titles = {'graphst': 'GraphST', 'stagate': 'STAGATE', 'spagcn': 'SpaGCN', 'bayesspace':'BayesSpace'}  # Pretty titles
    
    combined_data = []
    for method in methods:
        suffix = method_suffix[method]
        ari_scores = extract_ari_scores_all(technology, dataset, method, slice_num_list, suffix, hvg)
        # Convert the ARI scores to a DataFrame
        ari_df = pd.DataFrame.from_dict({(i, j): ari_scores[i][j]
                                         for i in ari_scores.keys()
                                         for j in ari_scores[i].keys()},
                                        orient='index')
    
        # Reset index to make normalization and slice number columns
        ari_df_reset = ari_df.reset_index()
        ari_df_reset.columns = ['Normalization', 'Slice_Num'] + [f"Run_{run}" for run in range(1, 6)]
    
        # Calculate the median ARI score across runs for each normalization and slice
        ari_df_reset['Median_ARI'] = ari_df_reset[[f"Run_{run}" for run in range(1, 6)]].median(axis=1)
    
        # Combine the original methods (stagate, spagcn, graphst) under "Original Method"
        ari_df_reset['Normalization'] = ari_df_reset['Normalization'].replace({
            'original_stagate': 'Original Method',
            'original_spagcn': 'Original Method',
            'original_graphst': 'Original Method',
            'original_bayesspace': 'Original Method'
        })
    
        # Replace other normalization names with friendly names
        name_dict = {
            "raw": "y",
            "raw_size": "y/s",
            "cpm": "CPM",
            "shifted_log": "log(y/s + 1)",
            "cpm_shifted_log": "log(CPM + 1)",
            "shifted_log_size": "log(y/s + 1)/u",
            "acosh": "acosh(2αy/s + 1)",
            "pseudo_shifted_log": "log(y/s + 1/4α)",
            "analytic_pearson_residual_clip": "Analytic Pearson (Clip)",
            "analytic_pearson_residual_noclip": "Analytic Pearson (No Clip)",
            "deseq2_log": "DESeq2 (log)",
            "deseq2_log1p": "DESeq2 (log1p)",
            "dino": "Dino",
            "scanpy_zheng": "scanpy Zheng",
            "scanpy_weinreb": "scanpy Weinreb",
            "scanpy_seurat": "scanpy Seurat",
            "scanpy_pearson_residual": "scanpy Pearson Residual",
            "tmm": "TMM",
            "normalisr": "Normalisr",
            "psinorm": "PsiNorm",
            "spanorm":"SpaNorm"
        }
        
        ari_df_reset['Normalization'] = ari_df_reset['Normalization'].replace(name_dict)
    
        # Add a column for the method to identify the method in the plot
        ari_df_reset['Method'] = method_titles[method]
        
        # Keep only necessary columns
        combined_data.append(ari_df_reset[['Normalization', 'Method', 'Median_ARI']])
    
    # Combine all the data into one DataFrame
    combined_df = pd.concat(combined_data)
    return combined_data, combined_df

# extract all the ari scores (including the raw counts and the original methods)
combined_data, combined_july = extract_combined('july')
df_aggregated = combined_july.groupby(["Method", "Normalization"], as_index=False)["Median_ARI"].median()
df_pivot = df_aggregated.pivot(index="Method", columns="Normalization", values="Median_ARI")
df_pivot = df_pivot.rename(columns={'Analytic Pearson Residual': 'Analytic Pearson'})
df_pivot.to_csv('/oscar/data/yma16/Project/spTransform/2.Downstream/results/summary/SDD/SDD_median_median_ARI_july_hvg3000.csv')

# aggregate the rank
df_pivot = df_pivot.drop(columns=["Original Method", "y",'Method']) # remove raw counts & original methods
df_pivot = df_pivot.apply(pd.to_numeric, errors='coerce')
ranked_df = df_pivot.rank(axis=1, method='average', ascending=True)
normalized_ranks = (ranked_df-1) / (len(ranked_df.columns)-1) # make the rank in range 0-1
total_rank = normalized_ranks.mean(axis=0)
total_rank_sorted = total_rank.sort_values(ascending=False)
total_rank_sorted.to_csv('/oscar/data/yma16/Project/spTransform/code/Manuscript/Version5/data/SDD_total_rank.csv')