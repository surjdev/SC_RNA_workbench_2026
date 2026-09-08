"""Reproducibility Acceptance Tests (NFR-1, Acceptance Criterion 3)

Verifies that running analysis functions twice with the same configuration and random seed
produces identical numerical and categorical results.
"""

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from workbench_utils import load_config
from workbench_utils.de import prepare_pydeseq2_data, run_pydeseq2
from workbench_utils.qc import calculate_qc_metrics, detect_doublets


@pytest.fixture
def reproducible_adata():
    rng = np.random.default_rng(42)
    n_cells, n_genes = 80, 100
    X = rng.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(np.float32)
    adata = sc.AnnData(
        X=X,
        obs=pd.DataFrame(
            {
                "condition": ["Control"] * 40 + ["Treated"] * 40,
                "batch": (["B1", "B2"] * 40),
            },
            index=[f"cell_{i}" for i in range(n_cells)],
        ),
        var=pd.DataFrame(
            {"gene_name": [f"gene_{i}" for i in range(n_genes)]},
            index=[f"gene_{i}" for i in range(n_genes)],
        ),
    )
    # Add mitochondrial genes
    adata.var_names = [f"MT-{i}" if i < 10 else f"GENE-{i}" for i in range(n_genes)]
    return adata


def test_qc_and_doublet_reproducibility(reproducible_adata):
    cfg = load_config()
    seed = cfg["project"]["random_seed"]

    # Run 1
    adata1 = reproducible_adata.copy()
    calculate_qc_metrics(adata1)
    adata1 = detect_doublets(adata1, random_state=seed, n_prin_comps=5)

    # Run 2
    adata2 = reproducible_adata.copy()
    calculate_qc_metrics(adata2)
    adata2 = detect_doublets(adata2, random_state=seed, n_prin_comps=5)

    # Assert identical results
    np.testing.assert_allclose(
        adata1.obs["doublet_score"].values,
        adata2.obs["doublet_score"].values,
        err_msg="Scrublet doublet scores differ between runs with same seed!",
    )
    pd.testing.assert_series_equal(
        adata1.obs["predicted_doublet"],
        adata2.obs["predicted_doublet"],
        check_names=True,
    )


def test_clustering_reproducibility(reproducible_adata):
    cfg = load_config()
    seed = cfg["project"]["random_seed"]

    def run_pipeline(ad):
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)
        sc.pp.highly_variable_genes(ad, n_top_genes=50)
        sc.pp.pca(ad, n_comps=10, random_state=seed)
        sc.pp.neighbors(ad, n_pcs=10, random_state=seed)
        sc.tl.leiden(ad, resolution=0.5, random_state=seed)
        return ad

    adata1 = run_pipeline(reproducible_adata.copy())
    adata2 = run_pipeline(reproducible_adata.copy())

    # PCA coordinates must match exactly
    np.testing.assert_allclose(
        adata1.obsm["X_pca"],
        adata2.obsm["X_pca"],
        atol=1e-6,
        err_msg="PCA coordinates differ between runs!",
    )
    # Leiden clusters must match exactly
    pd.testing.assert_series_equal(
        adata1.obs["leiden"],
        adata2.obs["leiden"],
        check_names=True,
    )


def test_pydeseq2_reproducibility(reproducible_adata):
    counts, meta = prepare_pydeseq2_data(
        reproducible_adata,
        design_factor="condition",
        layer=None,
        min_cells_per_gene=3,
    )

    stat_res1 = run_pydeseq2(
        counts_df=counts,
        clinical_df=meta,
        design_factors=["condition"],
        contrast=("condition", "Treated", "Control"),
        quiet=True,
    )

    stat_res2 = run_pydeseq2(
        counts_df=counts,
        clinical_df=meta,
        design_factors=["condition"],
        contrast=("condition", "Treated", "Control"),
        quiet=True,
    )

    # Numerical results must be identical
    np.testing.assert_allclose(
        stat_res1["log2FoldChange"].values,
        stat_res2["log2FoldChange"].values,
        rtol=1e-5,
        err_msg="PyDESeq2 log2FoldChange differs between identical runs!",
    )
