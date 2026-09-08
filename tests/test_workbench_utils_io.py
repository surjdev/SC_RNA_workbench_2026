import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from workbench_utils import export_for_seurat, load_h5ad, load_upstream_matrix, save_h5ad


@pytest.fixture
def sample_adata():
    np.random.seed(42)
    n_cells = 25
    n_genes = 50
    genes = [f"GENE_{i}" for i in range(40)] + [f"MT-{i}" for i in range(10)]
    cells = [f"cell_{i}" for i in range(n_cells)]
    X = sp.csr_matrix(np.random.poisson(lam=10, size=(n_cells, n_genes)).astype(np.float32))
    obs = pd.DataFrame(index=cells)
    obs["batch"] = ["batch_1"] * 12 + ["batch_2"] * 13
    var = pd.DataFrame(index=genes)
    var["gene_name"] = genes
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()
    adata.obsm["X_umap"] = np.random.randn(n_cells, 2)
    return adata


def test_save_and_load_h5ad_backed(sample_adata, tmp_path):
    out_file = tmp_path / "test_data.h5ad"
    save_h5ad(sample_adata, out_file)
    assert out_file.exists()

    # Test in-memory load
    adata_loaded = load_h5ad(out_file)
    assert adata_loaded.n_obs == sample_adata.n_obs
    assert adata_loaded.n_vars == sample_adata.n_vars

    # Test backed mode (NFR-3)
    adata_backed = load_h5ad(out_file, backed="r")
    assert adata_backed.isbacked
    assert adata_backed.n_obs == sample_adata.n_obs


def test_load_upstream_matrix(tmp_path):
    matrix_file = tmp_path / "counts.tsv"
    genes = [f"Gene_{i}" for i in range(10)]
    cells = [f"Cell_{j}" for j in range(5)]
    df = pd.DataFrame(np.random.poisson(lam=5, size=(10, 5)), index=genes, columns=cells)
    df.to_csv(matrix_file, sep="\t")

    adata = load_upstream_matrix(matrix_file, transpose=True, sparse=True)
    assert adata.n_obs == 5  # Cells
    assert adata.n_vars == 10  # Genes
    assert "counts" in adata.layers
    assert sp.issparse(adata.X)


def test_export_for_seurat(sample_adata, tmp_path):
    export_dir = tmp_path / "seurat_export"
    export_for_seurat(sample_adata, export_dir)

    assert (export_dir / "counts.csv").exists()
    assert (export_dir / "metadata.csv").exists()
    assert (export_dir / "embeddings.csv").exists()
