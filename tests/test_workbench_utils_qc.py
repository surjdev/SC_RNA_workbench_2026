import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from workbench_utils import (
    calculate_qc_metrics,
    detect_doublets,
    filter_cells,
    filter_genes,
    parse_plate_metadata,
)


@pytest.fixture
def mock_adata():
    np.random.seed(42)
    n_cells = 60
    n_genes = 85
    genes = (
        [f"GENE_{i}" for i in range(70)]
        + [f"MT-ND{i}" for i in range(10)]
        + [f"ERCC-{i:05d}" for i in range(5)]
    )
    cells = [f"Plate1_{chr(65 + (i // 12))}{(i % 12 + 1):02d}" for i in range(n_cells)]
    X = sp.csr_matrix(np.random.poisson(lam=15, size=(n_cells, n_genes)).astype(np.float32))
    obs = pd.DataFrame(index=cells)
    var = pd.DataFrame(index=genes)
    var["gene_name"] = genes
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()
    return adata


def test_qc_metrics_and_filtering(mock_adata):
    adata = calculate_qc_metrics(mock_adata)
    assert "n_genes_by_counts" in adata.obs
    assert "total_counts" in adata.obs
    assert "pct_counts_mito" in adata.obs
    assert "pct_counts_ercc" in adata.obs
    assert adata.var["mt"].sum() == 10
    assert adata.var["ercc"].sum() == 5

    # Filter cells with ERCC threshold
    filtered_cells = filter_cells(
        adata, min_genes=5, min_counts=20, max_pct_mito=50.0, max_pct_ercc=50.0
    )
    assert filtered_cells.n_obs > 0
    assert filtered_cells.n_obs <= adata.n_obs

    # Filter genes
    filtered_genes = filter_genes(filtered_cells, min_cells=2)
    assert filtered_genes.n_vars > 0


def test_parse_plate_metadata(mock_adata):
    adata = parse_plate_metadata(mock_adata)
    assert "plate" in adata.obs
    assert "well" in adata.obs
    assert "well_row" in adata.obs
    assert "well_col" in adata.obs
    assert adata.obs["plate"].iloc[0] == "Plate1"


def test_doublet_detection(mock_adata):
    adata = detect_doublets(mock_adata, expected_doublet_rate=0.05, random_state=42, n_prin_comps=5)
    assert "doublet_score" in adata.obs
    assert "predicted_doublet" in adata.obs
    assert len(adata.obs["doublet_score"]) == mock_adata.n_obs
