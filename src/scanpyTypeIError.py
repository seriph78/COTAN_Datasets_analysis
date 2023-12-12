#scanpyTypeIError.py

import sys
import numpy as np
import pandas as pd
import scanpy as sc
import os
import matplotlib.pyplot as plt
from scipy.stats import rankdata

def scanpyTypeIError(exprs, meta, feature_meta, embedding, mtCode, outDir):
  #exprs = sys.argv[1]
  #meta = sys.argv[2]
  #feature_meta = sys.argv[3]
  #embedding = sys.argv[4]
  
  #adata_seurat = sc.AnnData(X = r.exprs.T, obs = r.meta, var = r.feature_meta)
  adata_seurat = sc.AnnData(X = exprs.T, obs = meta, var = feature_meta)
  
  #adata_seurat.obsm['umap'] = r.embedding
  adata_seurat.obsm['umap'] = embedding
  
  adata_seurat.var_names_make_unique()
  adata_seurat
  
  matrix = adata_seurat.X
  matrix = matrix.todense()
  neg_variances = np.sort(-np.var(matrix, axis=0))
  sorted_log_variances = [np.log(-i) for i in neg_variances.T][:2000]
  #plt.scatter([i for i in range(len(sorted_log_variances))], sorted_log_variances)
  #plt.show()
  
  
  sc.pl.highest_expr_genes(adata_seurat, n_top=20, )
  
  
  sc.pp.filter_cells(adata_seurat, min_genes=200)
  sc.pp.filter_genes(adata_seurat, min_cells=3)
  
  adata_seurat.var['mt'] = adata_seurat.var_names.str.startswith(mtCode)  # annotate the group of mitochondrial genes as 'mt'
  sc.pp.calculate_qc_metrics(adata_seurat, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
  
  sc.pl.violin(adata_seurat, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],jitter=0.4, multi_panel=True)
  
  sc.pp.normalize_total(adata_seurat, target_sum=1e4)
  sc.pp.log1p(adata_seurat)
  
  sc.pp.highly_variable_genes(adata_seurat,
      n_top_genes=500, # chosen based on elbow plot above
      flavor='seurat_v3')
  
  sc.pl.highly_variable_genes(adata_seurat)
  
  adata_seurat.raw = adata_seurat
  adata_seurat = adata_seurat[:, adata_seurat.var.highly_variable]
  # sc.pp.regress_out(adata_seurat, ['total_counts', 'pct_counts_mt']) # analogous to Seurat's regressOut, not done in Seurat pipeline
  sc.pp.scale(adata_seurat, max_value=10) # FIXME: 10 in the tutorial, default None
  
  
  adata_seurat.obs['TestCl'] = adata_seurat.obs['TestCl'].astype('category')
  # Now run the differential expression analysis
  sc.tl.rank_genes_groups(adata_seurat, 'TestCl', method='wilcoxon')
  
  results = adata_seurat.uns['rank_genes_groups']
  groups = results['names'].dtype.names
  
  # Creating a DataFrame to hold results
  de_genes = pd.DataFrame({group: results['names'][group] for group in groups})
  
  # Adding adjusted p-values to the DataFrame
  for group in groups:
      de_genes[f'pvals_adj_{group}'] = results['pvals_adj'][group]
  
  # Filter genes with adjusted p-value < 0.05
  for group in groups:
      de_genes_filtered = de_genes[de_genes[f'pvals_adj_{group}'] < 0.05]
  
  # Display the filtered results
  print(de_genes_filtered)
  
  
  print(de_genes_filtered.size)
  de_genes_filtered.to_csv(outDir+"de_genes_filtered.csv")


#import scvi
