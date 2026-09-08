"""
Publication-Grade Visualization Module for Single-Cell Analysis.
Provides tailored palettes (Nature, Cell), violin QC plots, volcano plots for DE results,
and embedding visualization helpers.
Complies with scRNAseq_Workbench_Requirements.md (FR-6, NFR-4).
"""

from typing import Dict, List, Optional, Tuple

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

NATURE_PALETTE: List[str] = [
    "#E64B35",
    "#4DBBD5",
    "#00A087",
    "#3C5488",
    "#F39B7F",
    "#8491B4",
    "#91D1C2",
    "#DC0000",
    "#7E6148",
    "#B09C85",
]

CELL_PALETTE: List[str] = [
    "#20854E",
    "#0072B5",
    "#BC3C29",
    "#EEA236",
    "#6F99AD",
    "#FFDC00",
    "#79AF97",
    "#3B4992",
    "#E18727",
    "#20854E",
]


def set_publication_style(palette: str = "nature", dpi: int = 150) -> None:
    """Configure aesthetic publication styling for matplotlib and Scanpy figures."""
    sc.set_figure_params(dpi=dpi, frameon=False, vector_friendly=True, facecolor="white")
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Helvetica", "Arial"]
    plt.rcParams["axes.edgecolor"] = "#2b2b2b"
    plt.rcParams["axes.linewidth"] = 0.8
    plt.rcParams["xtick.major.width"] = 0.8
    plt.rcParams["ytick.major.width"] = 0.8


def plot_qc_violins(
    adata: ad.AnnData,
    keys: Tuple[str, ...] = ("n_genes_by_counts", "total_counts", "pct_counts_mito"),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Plot publication-quality violin distributions for single-cell QC metrics."""
    set_publication_style()

    available_keys = [k for k in keys if k in adata.obs]
    if not available_keys:
        raise KeyError(f"None of {keys} found in adata.obs!")

    fig, axes = plt.subplots(1, len(available_keys), figsize=(4 * len(available_keys), 4))
    if len(available_keys) == 1:
        axes = [axes]

    titles: Dict[str, str] = {
        "n_genes_by_counts": "Genes per Cell",
        "total_counts": "Total Reads / UMI",
        "pct_counts_mito": "Mitochondrial Reads (%)",
        "pct_counts_ribo": "Ribosomal Reads (%)",
        "doublet_score": "Doublet Score",
    }

    for i, key in enumerate(available_keys):
        sns.violinplot(
            y=adata.obs[key],
            ax=axes[i],
            color=NATURE_PALETTE[i % len(NATURE_PALETTE)],
            inner="quartile",
            linewidth=1.0,
        )
        axes[i].set_title(titles.get(key, key), fontsize=12, fontweight="bold")
        axes[i].set_ylabel("")

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight", dpi=300)
    return fig


def plot_volcano(
    de_df: pd.DataFrame,
    pval_col: str = "padj",
    fc_col: str = "log2FoldChange",
    gene_col: str = "gene",
    pval_threshold: float = 0.05,
    fc_threshold: float = 1.0,
    top_n_labels: int = 10,
    title: str = "Differential Expression Volcano Plot",
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Generate publication-ready Volcano Plot for PyDESeq2 / DE results.

    Points are colored by significance and fold-change threshold.
    """
    set_publication_style()
    df = de_df.copy()

    # Drop NA pvalues
    df = df.dropna(subset=[pval_col, fc_col])
    df["-log10(padj)"] = -np.log10(df[pval_col].replace(0, 1e-300))

    # Classify points
    conditions = [
        (df[pval_col] < pval_threshold) & (df[fc_col] >= fc_threshold),
        (df[pval_col] < pval_threshold) & (df[fc_col] <= -fc_threshold),
    ]
    choices = ["Up-regulated", "Down-regulated"]
    df["status"] = np.select(conditions, choices, default="Not Significant")

    color_map = {
        "Up-regulated": "#E64B35",
        "Down-regulated": "#4DBBD5",
        "Not Significant": "#B09C85",
    }

    fig, ax = plt.subplots(figsize=(7, 6))
    for status, grp in df.groupby("status"):
        ax.scatter(
            grp[fc_col],
            grp["-log10(padj)"],
            c=color_map.get(status, "#B09C85"),
            label=f"{status} ({len(grp)})",
            alpha=0.7,
            edgecolors="none",
            s=25,
        )

    # Threshold guidelines
    ax.axhline(-np.log10(pval_threshold), color="grey", linestyle="--", linewidth=0.8)
    ax.axvline(fc_threshold, color="grey", linestyle="--", linewidth=0.8)
    ax.axvline(-fc_threshold, color="grey", linestyle="--", linewidth=0.8)

    # Annotate top genes
    if gene_col in df.columns and top_n_labels > 0:
        top_genes = (
            df[df["status"] != "Not Significant"].sort_values(by=pval_col).head(top_n_labels)
        )
        for _, row in top_genes.iterrows():
            ax.annotate(
                row[gene_col],
                (row[fc_col], row["-log10(padj)"]),
                textcoords="offset points",
                xytext=(4, 4),
                fontsize=8,
                fontweight="bold",
            )

    ax.set_xlabel("log2 Fold Change", fontsize=11, fontweight="bold")
    ax.set_ylabel(f"-log10({pval_col})", fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.legend(frameon=False)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight", dpi=300)
    return fig
