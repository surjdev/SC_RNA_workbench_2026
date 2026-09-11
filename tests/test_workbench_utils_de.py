import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from workbench_utils import (
    RDependencyError,
    check_r_dependencies,
    deseq2_python,
    deseq2_r,
    prepare_pseudobulk_data,
    prepare_pydeseq2_data,
    run_de,
    run_pydeseq2,
)


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


def test_prepare_pseudobulk_data(de_adata):
    """Test data extraction and pseudobulk aggregation."""
    counts_df, clinical_df = prepare_pseudobulk_data(de_adata, design_factor="condition")
    assert counts_df.shape[0] == 30
    assert counts_df.shape[1] <= 20
    assert "condition" in clinical_df.columns

    # Test pseudo-bulk aggregation
    pb_counts, pb_meta = prepare_pseudobulk_data(
        de_adata, design_factor="condition", sample_key="sample_id"
    )
    assert pb_counts.shape[0] == 6
    assert "sample_id" in pb_meta.columns
    # Verify backward compatibility
    assert prepare_pydeseq2_data is prepare_pseudobulk_data


def test_deseq2_python_execution(de_adata):
    """Test pure-Python PyDESeq2 engine execution (FR-11)."""
    counts_df, clinical_df = prepare_pseudobulk_data(de_adata, design_factor="condition")
    results = deseq2_python(
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
    # Verify alias
    assert run_pydeseq2 is deseq2_python


def test_check_r_dependencies():
    """Test R dependency check function (FR-12, FR-13)."""
    is_avail, msg = check_r_dependencies(["DESeq2"])
    assert isinstance(is_avail, bool)
    if not is_avail:
        assert isinstance(msg, str)
        assert "install" in msg.lower() or "r-stats" in msg.lower()


def test_r_dependency_error_handling(de_adata):
    """Test that missing R/Bioconductor dependencies raise actionable RDependencyError (FR-13)."""
    counts_df, clinical_df = prepare_pseudobulk_data(de_adata, design_factor="condition")

    is_avail, _ = check_r_dependencies(["DESeq2"])
    if not is_avail:
        with pytest.raises(RDependencyError) as exc_info:
            deseq2_r(counts_df, clinical_df, design_factors="condition")
        err_msg = str(exc_info.value)
        assert "install" in err_msg.lower() or "r-stats" in err_msg.lower()
        # Ensure no raw cryptic traceback, but an informative message
        assert len(err_msg) > 20


def test_run_de_dispatcher_python(de_adata):
    """Test config-driven dispatcher for Python engine."""
    counts_df, clinical_df = prepare_pseudobulk_data(de_adata, design_factor="condition")
    results = run_de(
        counts_df=counts_df,
        clinical_df=clinical_df,
        design_factors="condition",
        contrast=("condition", "GroupA", "GroupB"),
        engine="pydeseq2",
        quiet=True,
    )
    assert not results.empty
    assert "gene" in results.columns


def test_run_de_dispatcher_invalid_engine(de_adata):
    """Test dispatcher rejects invalid engine with helpful ValueError."""
    counts_df, clinical_df = prepare_pseudobulk_data(de_adata, design_factor="condition")
    with pytest.raises(ValueError) as exc_info:
        run_de(
            counts_df=counts_df,
            clinical_df=clinical_df,
            design_factors="condition",
            engine="unknown_engine",
        )
    assert "Unsupported differential expression engine" in str(exc_info.value)


def test_deseq2_r_schema_and_mock(de_adata, monkeypatch):
    """Test deseq2_r converts R results into the exact expected pandas DataFrame schema."""
    import unittest.mock as mock

    counts_df, clinical_df = prepare_pseudobulk_data(de_adata, design_factor="condition")

    # Mock check_r_dependencies to return True
    monkeypatch.setattr("workbench_utils.de.check_r_dependencies", lambda pkgs: (True, None))

    fake_r_df = pd.DataFrame(
        {
            "baseMean": [10.5, 25.1],
            "log2FoldChange": [1.8, -2.1],
            "lfcSE": [0.3, 0.4],
            "stat": [6.0, -5.25],
            "pvalue": [1e-5, 2e-4],
            "padj": [1e-4, 1e-3],
            "gene": ["Gene_0", "Gene_1"],
        }
    )

    mock_robjects = mock.MagicMock()
    mock_robjects.r.return_value = "mock_r_res"
    mock_robjects.conversion.rpy2py.return_value = fake_r_df.copy()
    mock_robjects.conversion.py2rpy.return_value = "mock_rpy"
    mock_pandas2ri = mock.MagicMock()

    class MockLocalConverter:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_val, exc_tb):
            pass

    mock_converter = MockLocalConverter()

    mock_conversion = mock.MagicMock(
        localconverter=lambda *args: mock_converter,
        py2rpy=lambda x: x,
        rpy2py=lambda x: fake_r_df.copy(),
    )

    with mock.patch.dict(
        "sys.modules",
        {
            "rpy2": mock.MagicMock(),
            "rpy2.robjects": mock_robjects,
            "rpy2.robjects.pandas2ri": mock_pandas2ri,
            "rpy2.robjects.conversion": mock_conversion,
        },
    ):
        res = deseq2_r(
            counts_df,
            clinical_df,
            design_factors="condition",
            contrast=("condition", "GroupA", "GroupB"),
        )
        assert isinstance(res, pd.DataFrame)
        assert list(res["gene"]) == ["Gene_0", "Gene_1"]
        assert "padj" in res.columns
        assert "log2FoldChange" in res.columns
        assert "baseMean" in res.columns
