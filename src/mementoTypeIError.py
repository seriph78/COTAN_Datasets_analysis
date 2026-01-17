
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
from contextlib import nullcontext
from joblib import parallel_backend


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


def mementoDEA(dataRaw, genes,barcodes,code,
                   labels_df,
                   dirOut,                    
                   num_cpus,
                   num_boot,
                   capture_rate,
                   obs_key="fake_cluster"
                   ):
    
    #dataRaw = mmread(inDir+name).tocsc()
    num_cpus = int(num_cpus)
    num_boot = int(num_boot)
    capture_rate = float(capture_rate)
    
    print(dataRaw.shape, dataRaw.nnz, type(dataRaw))
    print("New dea new")

    os.makedirs(dirOut, exist_ok=True)
    X = dataRaw.T.tocsr() if sparse.issparse(dataRaw) else np.asarray(dataRaw).T
    adata = ad.AnnData(X=X)
    adata.var_names = pd.Index([str(g) for g in genes], dtype=str)
    adata.obs_names = pd.Index([str(b) for b in barcodes], dtype=str)


    #adata =add_fake_cluster_by_library_size(adata, perc = percentage)

    # ---- align labels safely by barcode
    if not {"barcode", obs_key}.issubset(labels_df.columns):
        raise ValueError(f"labels_df must have columns: 'barcode' and '{obs_key}'")

    label_map = (
        labels_df[["barcode", obs_key]]
        .astype(str)
        .drop_duplicates(subset="barcode")
        .set_index("barcode")[obs_key]
    )

    adata.obs[obs_key] = label_map.reindex(adata.obs_names).values
    if pd.isna(adata.obs[obs_key]).any():
        missing = int(pd.isna(adata.obs[obs_key]).sum())
        raise ValueError(f"{missing} barcodes in adata are missing in labels_df (mismatch/order problem).")

    adata.obs[obs_key] = adata.obs[obs_key].astype("category")
    clusters = adata.obs[obs_key].cat.categories.tolist()

    all_1d_results = []
    for cl in clusters:
        col = f"is_cluster_{cl}"
        adata.obs[col] = (adata.obs[obs_key].astype(str) == str(cl)).astype(int)

        # Force threads instead of loky processes (MUCH safer with reticulate on Windows)
        with parallel_backend("threading", n_jobs=num_cpus):
            # res = memento.binary_test_1d(
            #     adata=adata,
            #     capture_rate=capture_rate,
            #     treatment_col=col,
            #     num_cpus=num_cpus,
            #     num_boot=num_boot,
            # )
            res = binary_test_1d_all_genes(
                adata=adata,
                capture_rate=capture_rate,
                treatment_col=col,
                num_cpus=num_cpus,
                num_boot=num_boot,
                filter_mean_thresh=0.0,
                min_perc_group=0.0,   # or 0.0 if you only want genes expressed in ≥1 group
            )
        res["cluster"] = cl
        all_1d_results.append(res)

    df_merged = pd.concat(all_1d_results, axis=0, ignore_index=True)

    # Bonferroni over genes
    m = adata.n_vars
    df_merged["de_pval_adj"] = (df_merged["de_pval"].astype(float) * m).clip(upper=1.0)

    out_path = os.path.join(dirOut, f"{code}_Memento_DEA_genes.csv")
    df_merged.to_csv(out_path, index=False)

    return df_merged


def binary_test_1d_all_genes(
    adata,
    capture_rate,
    treatment_col,
    num_cpus,
    num_boot=5000,
    verbose=1,
    replicates=(),
    filter_mean_thresh=0.0,
    min_perc_group=0.0,     # IMPORTANT: ≥1 group expressed
    min_cell_count=1,       # IMPORTANT: don't drop small groups
):
    adata = adata.copy()
    original_genes = adata.var_names.astype(str).to_list()

    adata.obs["capture_rate"] = float(capture_rate)

    memento.setup_memento(
        adata,
        q_column="capture_rate",
        filter_mean_thresh=filter_mean_thresh,
        min_cell_count=min_cell_count,
    )  # setup_memento defaults include filter_mean_thresh=0.07 

    memento.create_groups(adata, label_columns=[treatment_col] + list(replicates))

    # IMPORTANT: allow subsetting to "testable" genes so ht_1d doesn't crash
    memento.compute_1d_moments(
        adata,
        min_perc_group=min_perc_group,
        filter_genes=True,
    )  # overall_gene_mask uses gene_filter_rate > min_perc_group 

    sample_meta = memento.get_groups(adata)[[treatment_col]]

    memento.ht_1d_moments(
        adata,
        treatment=sample_meta,
        num_boot=num_boot,
        verbose=verbose,
        num_cpus=num_cpus,
    )

    res = memento.get_1d_ht_result(adata)

    # Expand back to all genes (for ROC)
    full = pd.DataFrame({"gene": original_genes}).merge(res, on="gene", how="left")

    # For ROC: untested genes = non-sig
    for c in ["de_pval", "dv_pval"]:
        if c in full.columns:
            full[c] = pd.to_numeric(full[c], errors="coerce").fillna(1.0)

    return full