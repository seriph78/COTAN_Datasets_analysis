#scanpyFDR.py
# This script is usabe only for two clusters situation!

import sys
import numpy as np
import pandas as pd
import scanpy as sc
import os
import matplotlib.pyplot as plt
import scipy

from scipy.stats import rankdata
from matplotlib.backends.backend_pdf import PdfPages

def scanpyFDR(exprs, meta, feature_meta,  mtCode, outDir, sample):
  #exprs = sys.argv[1]
  #meta = sys.argv[2]
  #feature_meta = sys.argv[3]
  #embedding = sys.argv[4]
  
  #adata_seurat = sc.AnnData(X = r.exprs.T, obs = r.meta, var = r.feature_meta)
  adata = sc.AnnData(X = exprs.T, obs = meta, var = feature_meta)
  
  #adata_seurat.obsm['umap'] = r.embedding
  #adata_seurat.obsm['umap'] = embedding
  
  adata.var_names_make_unique()
  #adata_seurat
  
  matrix = adata.X
  matrix = matrix.todense()
  neg_variances = np.sort(-np.var(matrix, axis=0))
  sorted_log_variances = [np.log(-i) for i in neg_variances.T][:2000]
  plt.scatter([i for i in range(len(sorted_log_variances))], sorted_log_variances)
  scatter = plt.show()
  
  pp = PdfPages(outDir+sample+'_plots.pdf')
  pp.savefig(scatter)
  
  pl_topGenes = sc.pl.highest_expr_genes(adata, n_top=20, )
  pp.savefig(pl_topGenes)
  
  sc.pp.filter_cells(adata, min_genes=20)
  sc.pp.filter_genes(adata, min_cells=3)
  
  adata.var['mt'] = adata.var_names.str.startswith(mtCode)  # annotate the group of mitochondrial genes as 'mt'
  sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
  
  violin = sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],jitter=0.4, multi_panel=True)
  pp.savefig(violin)
  
  pl1 = sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
  pl2 = sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
  pp.savefig(pl1)
  pp.savefig(pl2)
  
  sc.pp.normalize_total(adata, target_sum=1e4)
  sc.pp.log1p(adata)
  
  sc.pp.highly_variable_genes(adata,min_mean=0.0125, max_mean=3, min_disp=0.5)
      #n_top_genes=500)#, # chosen based on elbow plot above
 #     flavor='seurat_v3')
  
  sc.pl.highly_variable_genes(adata)
  hi_var_genes = plt.show()
  pp.savefig(hi_var_genes)
  
  adata.raw = adata
  #print(adata_seurat)
  adata = adata[:, adata.var.highly_variable]
  #print(adata_seurat)
  sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt']) # analogous to Seurat's regressOut, not done in Seurat pipeline
  sc.pp.scale(adata, max_value=10) # FIXME: 10 in the tutorial, default None
  
  
  adata.obs['TestCl'] = adata.obs['TestCl'].astype('category')
  
  # Now run the differential expression analysis
  sc.tl.rank_genes_groups(adata, 'TestCl', method='wilcoxon')
  
  # Extracting the differential expression results
  results = adata.uns['rank_genes_groups']
  #print(results)
  #print(results_pre)
  groups = results['names'].dtype.names
  #print(groups[0])

  # Create a DataFrame for the log fold changes, p-values, and adjusted p-values
  genes = results['names'][groups[0]]  # Adjust groups[0] to your cluster name
  logfoldchanges = results['logfoldchanges'][groups[0]]
  score = results['scores'][groups[0]]
  pvals = results['pvals'][groups[0]]
  pvals_adj = results['pvals_adj'][groups[0]]
  clusters = np.repeat("cl" + groups[0], len(genes))

  for cl in range(1, len(groups)):
      genes = np.append(genes, results['names'][groups[cl]])  # Adjust groups[cl] to your cluster name
      logfoldchanges = np.append(logfoldchanges, results['logfoldchanges'][groups[cl]])
      score = np.append(score, results['scores'][groups[cl]])
      pvals = np.append(pvals, results['pvals'][groups[cl]])
      pvals_adj = np.append(pvals_adj, results['pvals_adj'][groups[cl]])
      clusters = np.append(clusters, np.repeat("cl" + groups[cl], len(results['names'][groups[cl]])))

  df = pd.DataFrame({
      'genes': genes,
      'log_fold_change': logfoldchanges,
      'score': score,
      'pval': pvals,
      'pval_adj': pvals_adj,
      'clusters': clusters
  })
    
  # Calculating average expression for each cluster
  # Replace groups[0] and groups[1] with your actual cluster names
  
  # avg_exp_cl1 = np.mean(adata[adata.obs['TestCl'] == int(float(groups[0])), :].X, axis=0)
  # avg_exp_cl2 = np.mean(adata[adata.obs['TestCl'] == int(float(groups[1])), :].X, axis=0)
  # 
  # #print(adata_seurat[adata_seurat.obs['TestCl'] == int(float(groups[0])), :].X)
  #  #print(avg_exp_cl1.size())
  # #print(df.head())
  # # If the AnnData object's X is a sparse matrix, convert it to a dense format
  # # Convert the average expressions to dense format if they are sparse
  # 
  # avg_exp_cl1_dense = avg_exp_cl1.A1 if scipy.sparse.issparse(avg_exp_cl1) else avg_exp_cl1
  # avg_exp_cl2_dense = avg_exp_cl2.A1 if scipy.sparse.issparse(avg_exp_cl2) else avg_exp_cl2
  # 
  # # Create a DataFrame from the average expressions
  # avg_exp_df = pd.DataFrame({
  #     'avg_expr_cl1': avg_exp_cl1_dense, 
  #     'avg_expr_cl2': avg_exp_cl2_dense
  # }, index=adata_seurat.var_names)  
  # #print(avg_exp_df.shape)
  # #print(df.shape)
  # 
  # # Merge this DataFrame with df
  # df_merged = df.merge(avg_exp_df, left_on='genes', right_index=True)
  # print(df_merged.head())
  # # Set gene names as the index
  # df_merged.set_index('genes', inplace=True)
  # 
  #df_merged.to_csv(outDir+sample+"_Scampy_de_genes_"+str(percentage)+".csv")
  df.to_csv(outDir+sample+"_Scampy_DEA_all_genes.csv")

  #################
  # 
  # # Creating a DataFrame to hold results
  # de_genes = pd.DataFrame({group: results['names'][group] for group in groups})
  # 
  # # Adding adjusted p-values to the DataFrame
  # for group in groups:
  #     de_genes[f'pvals_adj_{group}'] = results['pvals_adj'][group]
  # 
  # # Filter genes with adjusted p-value < 0.05
  # # for group in groups:
  # #     de_genes_filtered = de_genes[de_genes[f'pvals_adj_{group}'] < 0.05]
  # # 
  # # Display the filtered results
  # #print(de_genes)
  # 
  # 
  # #print(de_genes.size)
  # de_genes.to_csv(outDir+sample+"de_genes_filtered.csv")

  pp.close()
#import scvi
