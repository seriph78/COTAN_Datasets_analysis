
import re
import os
import numpy as np
import pandas as pd
import anndata as ad
from io import StringIO
from scipy import sparse
from scipy.io import mmread
import anndata as ad
import memento

def add_fake_cluster_by_library_size(
    adata,
    perc: float,
    use_layer: str | None = None,
    obs_key: str = "fake_cluster",
):
    if not (0 < perc <= 1):
        raise ValueError("perc must be in (0, 1].")

    X = adata.layers[use_layer] if use_layer is not None else adata.X

    lib = np.asarray(X.sum(axis=1)).ravel() if sparse.issparse(X) else X.sum(axis=1)

    n = adata.n_obs
    k = int(np.rint(n * perc))
    k = max(1, min(k, n))

    top_idx = np.argsort(-lib)[:k]

    labels = np.full(n, 2, dtype=int)   # default cluster 2
    labels[top_idx] = 1                 # top perc -> cluster 1

    adata.obs[obs_key] = labels.astype(str)  # or .astype("category") with pandas
    adata.obs[obs_key] = adata.obs[obs_key].astype("category")

    return adata


def mementoTypeIError(name, 
                   inDir, 
                   dirOut, 
                   percentage,
                   num_cpus,
                   num_boot,
                   capture_rate
                   ):
    
    dataRaw = mmread(inDir+name).tocsc()
    print(dataRaw.shape, dataRaw.nnz, type(dataRaw))
    code = name.removesuffix('.mtx')

    genes = pd.read_csv(inDir+code+"genes.txt", header=None)[0].astype(str).to_list()
    barcodes = pd.read_csv(inDir+code+"barcodes.txt", header=None)[0].astype(str).to_list()

    adata = ad.AnnData(X=dataRaw.T.tocsr())
    adata.var_names = genes
    adata.obs_names = barcodes

    adata =add_fake_cluster_by_library_size(adata, perc = percentage)

    # capture_rate is the rough estimate of the overall UMI efficiency across both
    # sampling and sequencing. If s is the sequencing saturation, multiply s by
    # 0.07 for 10X v1, 0.15 for v2, and 0.25 for v3. This allows you to enter
    # different numbers for each batch, which likely have different saturation
    # numbers. This will NOT account for wildly different sequencing scenarios.
    capture_rate = capture_rate#0.25 # TODO: set properly
    num_cpus = num_cpus
    num_boot = num_boot #5000 # controls how many bootstrap iterations to perform, higher = more stable

    clusters = adata.obs["fake_cluster"].cat.categories.tolist()

    all_1d_results = []
    for cl in clusters:
        print("cluster "+cl)
        col = f"is_cluster_{cl}"
        adata.obs[col] = (adata.obs["fake_cluster"] == cl).astype(int)
        # run memento simple-binary DE (cluster vs rest)
        
        res = memento.binary_test_1d(
            adata=adata,
            capture_rate=capture_rate,
            treatment_col=col,
            num_cpus=num_cpus,
            num_boot=num_boot
        )
        res["cluster"] = cl
    #    res.to_csv(os.path.join(out_dir, f"memento_de_cluster_{cl}_vs_rest_1d.csv"), index=False)
        all_1d_results.append(res)
    df_merged = pd.concat(all_1d_results, axis=0, ignore_index=True)
    df_merged["de_pval_adj"] = np.minimum(df_merged["de_pval"].to_numpy(dtype=float) * adata.shape[1], 1.0)  # Bonferroni correction
    df_merged.to_csv(dirOut+code+"_memento"+"_"+str(percentage)+".csv", index=False)

    return df_merged