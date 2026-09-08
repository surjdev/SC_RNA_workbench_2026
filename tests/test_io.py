import anndata as ad
import numpy as np
import pandas as pd
import pytest

from sc_workbench.io import export_for_seurat, load_upstream_matrix, save_h5ad


@pytest.fixture
def sample_tsv(tmp_path):
    # Create mock Gene x Cell TSV
    genes = ["GAPDH", "ACTB", "TP53", "CD4"]
    cells = ["cell_1", "cell_2", "cell_3"]
    data = np.array([[100, 120, 80], [50, 60, 45], [10, 15, 8], [5, 8, 2]])
    df = pd.DataFrame(data, index=genes, columns=cells)
    tsv_p = tmp_path / "mock_counts.tsv"
    df.to_csv(tsv_p, sep="\t")
    return tsv_p


def test_load_upstream_matrix(sample_tsv):
    adata = load_upstream_matrix(sample_tsv, transpose=True)
    assert adata.n_obs == 3  # 3 cells
    assert adata.n_vars == 4  # 4 genes
    assert "counts" in adata.layers
    assert "gene_name" in adata.var


def test_save_h5ad(sample_tsv, tmp_path):
    adata = load_upstream_matrix(sample_tsv)
    out_p = tmp_path / "test.h5ad"
    save_h5ad(adata, out_p)
    assert out_p.exists()

    loaded = ad.read_h5ad(out_p)
    assert loaded.n_obs == adata.n_obs
    assert loaded.n_vars == adata.n_vars


def test_export_for_seurat(sample_tsv, tmp_path):
    adata = load_upstream_matrix(sample_tsv)
    adata.obsm["X_umap"] = np.random.randn(adata.n_obs, 2)
    export_dir = tmp_path / "seurat_out"
    export_for_seurat(adata, export_dir)
    assert (export_dir / "counts.csv").exists()
    assert (export_dir / "metadata.csv").exists()
    assert (export_dir / "embeddings_umap.csv").exists()
