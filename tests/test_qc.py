import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from sc_workbench.qc import calculate_qc_metrics, filter_cells, filter_genes


@pytest.fixture
def sample_adata():
    np.random.seed(42)
    n_cells = 20
    n_genes = 50
    # Include mitochondrial and ribosomal genes
    genes = (
        [f"GENE_{i}" for i in range(40)]
        + [f"MT-ATP{i}" for i in range(5)]
        + [f"RPS{i}" for i in range(5)]
    )
    cells = [f"cell_{i}" for i in range(n_cells)]
    X = sp.csr_matrix(np.random.poisson(lam=5, size=(n_cells, n_genes)).astype(np.float32))
    obs = pd.DataFrame(index=cells)
    var = pd.DataFrame(index=genes)
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()
    return adata


def test_calculate_qc_metrics(sample_adata):
    adata = calculate_qc_metrics(sample_adata)
    assert "total_counts" in adata.obs
    assert "n_genes_by_counts" in adata.obs
    assert "pct_counts_mito" in adata.obs
    assert "pct_counts_ribo" in adata.obs
    assert (adata.obs["total_counts"] > 0).all()


def test_filter_cells_and_genes(sample_adata):
    adata = calculate_qc_metrics(sample_adata)
    filtered = filter_cells(adata, min_genes=10, min_counts=50, max_pct_mito=30.0)
    assert filtered.n_obs <= adata.n_obs

    filtered_genes = filter_genes(filtered, min_cells=2)
    assert filtered_genes.n_vars <= adata.n_vars
