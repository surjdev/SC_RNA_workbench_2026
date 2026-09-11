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


def test_load_smartseq2_matrix(tmp_path):
    from workbench_utils import load_smartseq2_matrix

    counts_file = tmp_path / "smartseq2_counts.tsv"
    meta_file = tmp_path / "plate_meta.csv"

    genes = ["GAPDH", "ACTB", "ERCC-00001", "ERCC-00002"]
    cells = ["Plate1_A01", "Plate1_A02", "Plate2_B05"]

    counts_df = pd.DataFrame(
        np.array(
            [
                [1000, 1500, 2000],
                [800, 900, 1100],
                [50, 60, 40],
                [30, 40, 20],
            ]
        ),
        index=genes,
        columns=cells,
    )
    counts_df.to_csv(counts_file, sep="\t")

    meta_df = pd.DataFrame(
        {"treatment": ["Control", "Control", "Treated"], "donor": ["D1", "D1", "D2"]},
        index=cells,
    )
    meta_df.to_csv(meta_file)

    adata = load_smartseq2_matrix(counts_file, metadata_path=meta_file, transpose=True)
    assert adata.n_obs == 3
    assert adata.n_vars == 4
    assert adata.var["is_ercc"].sum() == 2
    assert "treatment" in adata.obs
    assert "plate" in adata.obs
    assert "well_row" in adata.obs
    assert adata.obs.loc["Plate1_A01", "well_row"] == "A"
    assert adata.obs.loc["Plate1_A01", "well_col"] == 1
