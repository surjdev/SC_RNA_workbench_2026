import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from workbench_utils import (
    plot_plate_layout,
    plot_qc_violins,
    plot_volcano,
    set_publication_style,
)


@pytest.fixture
def dummy_adata():
    obs = pd.DataFrame(
        {
            "n_genes_by_counts": [250, 400, 600],
            "total_counts": [1000, 2000, 3000],
            "pct_counts_mito": [2.5, 4.0, 5.5],
            "pct_counts_ercc": [1.0, 2.0, 3.0],
            "plate": ["Plate1", "Plate1", "Plate1"],
            "well_row": ["A", "B", "C"],
            "well_col": [1, 2, 3],
        },
        index=["Plate1_A01", "Plate1_B02", "Plate1_C03"],
    )
    return ad.AnnData(X=np.zeros((3, 3)), obs=obs)


def test_set_publication_style():
    set_publication_style(font_scale=1.0, style="whitegrid")
    assert plt.rcParams["axes.linewidth"] == 0.8


def test_plot_qc_violins(dummy_adata):
    fig = plot_qc_violins(dummy_adata)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_plate_layout(dummy_adata):
    fig = plot_plate_layout(dummy_adata, plate_id="Plate1", color_key="total_counts")
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_volcano():
    df = pd.DataFrame(
        {
            "gene": ["GeneA", "GeneB", "GeneC"],
            "log2FoldChange": [2.5, -3.0, 0.2],
            "padj": [0.001, 0.005, 0.8],
        }
    )
    fig = plot_volcano(df)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)
