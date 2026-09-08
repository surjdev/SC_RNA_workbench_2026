import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from sc_workbench import SingleCellWorkbench


@pytest.fixture
def test_wb():
    np.random.seed(42)
    n_cells = 30
    n_genes = 60
    genes = [f"GENE_{i}" for i in range(50)] + [f"MT-{i}" for i in range(10)]
    cells = [f"cell_{i}" for i in range(n_cells)]
    X = sp.csr_matrix(np.random.poisson(lam=10, size=(n_cells, n_genes)).astype(np.float32))
    obs = pd.DataFrame(index=cells)
    var = pd.DataFrame(index=genes)
    var["gene_name"] = genes
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()
    return SingleCellWorkbench(adata)


def test_full_fluent_pipeline(test_wb):
    wb = test_wb
    # 1. QC & Filter
    wb.calculate_qc().filter_cells(min_genes=5, min_counts=20, max_pct_mito=40.0)
    assert wb.n_cells > 0

    # 2. Preprocess
    wb.normalize(target_sum=1e4).select_hvg(n_top_genes=30).scale()
    assert "scaled" in wb.adata.layers

    # 3. Reduction & kNN
    wb.run_pca(n_comps=10).compute_neighbors(n_neighbors=5, n_pcs=10).run_umap()
    assert "X_pca" in wb.adata.obsm
    assert "X_umap" in wb.adata.obsm

    # 4. Clustering & Markers
    wb.cluster(resolution=0.5, method="leiden")
    assert "leiden" in wb.adata.obs

    markers_df = wb.find_markers(groupby="leiden", n_genes=5)
    assert not markers_df.empty
    assert "cluster" in markers_df.columns
    assert "gene" in markers_df.columns
