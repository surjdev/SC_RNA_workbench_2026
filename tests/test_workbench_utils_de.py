import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from workbench_utils import prepare_pydeseq2_data, run_pydeseq2


@pytest.fixture
def de_adata():
    np.random.seed(42)
    n_cells = 30
    n_genes = 20
    genes = [f"Gene_{i}" for i in range(n_genes)]
    cells = [f"cell_{i}" for i in range(n_cells)]

    # Make differential expression in first 5 genes
    X = np.random.poisson(lam=10, size=(n_cells, n_genes)).astype(np.float32)
    X[:15, :5] += 20  # Group A has higher expression

    obs = pd.DataFrame(index=cells)
    obs["condition"] = ["GroupA"] * 15 + ["GroupB"] * 15
    obs["sample_id"] = [f"sample_{i % 6}" for i in range(n_cells)]
    var = pd.DataFrame(index=genes)
    var["gene_name"] = genes

    adata = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)
    adata.layers["counts"] = sp.csr_matrix(X)
    return adata


def test_prepare_pydeseq2_data(de_adata):
    counts_df, clinical_df = prepare_pydeseq2_data(de_adata, design_factor="condition")
    assert counts_df.shape[0] == 30
    assert counts_df.shape[1] <= 20
    assert "condition" in clinical_df.columns

    # Test pseudo-bulk aggregation
    pb_counts, pb_meta = prepare_pydeseq2_data(
        de_adata, design_factor="condition", sample_key="sample_id"
    )
    assert pb_counts.shape[0] == 6
    assert "sample_id" in pb_meta.columns


def test_run_pydeseq2_execution(de_adata):
    counts_df, clinical_df = prepare_pydeseq2_data(de_adata, design_factor="condition")
    results = run_pydeseq2(
        counts_df=counts_df,
        clinical_df=clinical_df,
        design_factors="condition",
        contrast=("condition", "GroupA", "GroupB"),
        quiet=True,
    )
    assert not results.empty
    assert "gene" in results.columns
    assert "log2FoldChange" in results.columns
    assert "padj" in results.columns
