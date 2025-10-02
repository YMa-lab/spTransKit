#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import argparse
import os

def main():
    start_time = time.time()
    parser = argparse.ArgumentParser(description='downstream')
    parser.add_argument('--technology', type=str, default='10xVisium', help='techonology name')
    parser.add_argument('--dataset', type=str, default='DLPFC', help='dataset setting')
    parser.add_argument('--method', type=str, default='graphst')
    parser.add_argument('--graphst-radius', type=int, default=50)
    parser.add_argument('--stagate-cutoff', type=int, default=6)
    parser.add_argument('--hvg_type', type=str, default=None)
    args = parser.parse_args()
    
    
    file_dic = '/oscar/data/yma16/Project/spTransform/0.NormalizedCounts/' + args.technology + '/'
    file_dic_gt = '/oscar/data/yma16/Project/spTransform/2.Downstream/datasets/' + args.technology + '/' + args.dataset + '/'
    if args.dataset == 'DLPFC':
        slice_num_list = ['151508', '151509', '151507', '151671', '151676', '151670', '151669', '151510', '151675', '151672', '151673', '151674']
        
    if args.method == 'graphst':
        import torch
        from GraphST_pipeline import run_graphst, chech_fileexist
        filepath_output = f'/oscar/data/yma16/Project/spTransform/2.Downstream/results/{args.technology}/{args.dataset}_hvg_{args.hvg_type}_sdd/'
        
        args.cuda = torch.cuda.is_available()
        args.device = torch.device("cuda" if args.cuda else "cpu")
        for slice_num in slice_num_list:
            print('slice_num' + str(slice_num), flush=True)
            filepath_output_cp = filepath_output + slice_num + '/' + args.method + '/'
            ### current filepath_output is '/oscar/data/yma16/Project/spTransform/2.Downstream/results/10xVisium/DLPFC_hvg/slice_num/graphst/'
            if not os.path.exists(filepath_output_cp):
                os.makedirs(filepath_output_cp)
            if chech_fileexist(args, filepath_output_cp, slice_num) == False:
                run_graphst(args, file_dic, file_dic_gt, slice_num, filepath_output_cp, args.hvg_type, radius = args.graphst_radius)
    elif args.method == 'stagate':
        from STAGATE_pipeline import run_stagate, chech_fileexist_stagate
        filepath_output = f'/oscar/data/yma16/Project/spTransform/2.Downstream/results/{args.technology}/{args.dataset}_hvg_{args.hvg_type}_sdd/'
        for slice_num in slice_num_list:
            print('slice_num' + str(slice_num), flush=True)
            filepath_output_cp = filepath_output + slice_num + '/' + args.method + '/'
            if not os.path.exists(filepath_output_cp):
                os.makedirs(filepath_output_cp)
            if chech_fileexist_stagate(args, filepath_output_cp, slice_num) == False:
                run_stagate(args, file_dic, file_dic_gt, slice_num, filepath_output_cp, args.hvg_type, k_cutoff=args.stagate_cutoff)
    elif args.method == 'spagcn':
        from SpaGCN_pipeline import run_spagcn, chech_fileexist_spagcn
        filepath_output = f'/oscar/data/yma16/Project/spTransform/2.Downstream/results/{args.technology}/{args.dataset}_hvg_{args.hvg_type}_sdd/'
        for slice_num in slice_num_list:
            print('slice_num' + str(slice_num), flush=True)
            filepath_output_cp = filepath_output + slice_num + '/' + args.method + '/'
            if not os.path.exists(filepath_output_cp):
                os.makedirs(filepath_output_cp)
            if chech_fileexist_spagcn(args, filepath_output_cp, slice_num) == False:
                run_spagcn(args, file_dic, file_dic_gt, slice_num, filepath_output_cp, args.hvg_type)
                
        
    print("--- %s seconds ---" % (time.time() - start_time),flush=True)
    
    
if __name__ == '__main__':
    main()


