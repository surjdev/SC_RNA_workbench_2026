"""
Publication-Grade Visualization Module for Single-Cell Analysis.
Provides tailored palettes (Nature, Cell, Science), violin QC metrics, UMAPs,
volcano plots, and cluster dotplots.
"""

from typing import Dict, List, Optional, Tuple, Union

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

# Tailored color palettes
NATURE_PALETTE = [
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

CELL_PALETTE = [
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


def set_publication_style(palette: str = "nature", dpi: int = 150):
    """Set aesthetic publication styling for matplotlib and Scanpy."""
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

    fig, axes = plt.subplots(1, len(keys), figsize=(4 * len(keys), 4))
    if len(keys) == 1:
        axes = [axes]

    titles = {
        "n_genes_by_counts": "Genes per Cell",
        "total_counts": "Total Reads / UMI",
        "pct_counts_mito": "Mitochondrial Reads (%)",
        "doublet_score": "Doublet Score",
    }

    for i, key in enumerate(keys):
        if key in adata.obs:
            sns.violinplot(
                y=adata.obs[key],
                ax=axes[i],
                color=NATURE_PALETTE[i % len(NATURE_PALETTE)],
                inner="quartile",
                cut=0,
            )
            sns.stripplot(
                y=adata.obs[key], ax=axes[i], color="black", alpha=0.2, jitter=0.2, size=2
            )
            axes[i].set_title(titles.get(key, key), fontweight="bold", fontsize=11)
            axes[i].set_ylabel("")

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight", dpi=300)
    return fig


def plot_umap(
    adata: ad.AnnData,
    color: Union[str, List[str]] = "leiden",
    palette: Optional[List[str]] = None,
    ncols: int = 3,
    save_path: Optional[str] = None,
):
    """Plot 2D UMAP projection with custom styling."""
    set_publication_style()
    colors = [color] if isinstance(color, str) else color
    pal = palette or NATURE_PALETTE

    sc.pl.umap(adata, color=colors, palette=pal, ncols=ncols, show=False, frameon=False)
    if save_path:
        plt.savefig(save_path, bbox_inches="tight", dpi=300)


def plot_marker_dotplot(
    adata: ad.AnnData,
    markers: Union[Dict[str, List[str]], List[str]],
    groupby: str = "leiden",
    save_path: Optional[str] = None,
):
    """Generate a publication-grade DotPlot of cluster biomarkers."""
    set_publication_style()
    sc.pl.dotplot(
        adata, var_names=markers, groupby=groupby, standard_scale="var", cmap="Blues", show=False
    )
    if save_path:
        plt.savefig(save_path, bbox_inches="tight", dpi=300)


def plot_volcano(
    de_df: pd.DataFrame,
    title: str = "Volcano Plot",
    logfc_thresh: float = 1.0,
    pval_thresh: float = 0.05,
    top_n_labels: int = 10,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Generate a Volcano Plot highlighting significant differentially expressed genes."""
    set_publication_style()

    df = de_df.copy()
    df["neg_log10_pval"] = -np.log10(np.maximum(df["pvals_adj"], 1e-300))

    # Classify points
    conditions = [
        (df["logfoldchange"] >= logfc_thresh) & (df["pvals_adj"] <= pval_thresh),
        (df["logfoldchange"] <= -logfc_thresh) & (df["pvals_adj"] <= pval_thresh),
    ]
    df["significance"] = np.select(conditions, ["Up", "Down"], default="Not Significant")

    fig, ax = plt.subplots(figsize=(6, 5))
    palette_dict = {"Up": "#E64B35", "Down": "#4DBBD5", "Not Significant": "#D1D5DB"}

    sns.scatterplot(
        data=df,
        x="logfoldchange",
        y="neg_log10_pval",
        hue="significance",
        palette=palette_dict,
        alpha=0.8,
        s=25,
        ax=ax,
        edgecolor=None,
    )

    # Threshold lines
    ax.axvline(x=logfc_thresh, color="gray", linestyle="--", linewidth=0.8)
    ax.axvline(x=-logfc_thresh, color="gray", linestyle="--", linewidth=0.8)
    ax.axhline(y=-np.log10(pval_thresh), color="gray", linestyle="--", linewidth=0.8)

    # Label top genes
    top_genes = df[df["significance"] != "Not Significant"].nlargest(top_n_labels, "neg_log10_pval")
    for _, row in top_genes.iterrows():
        ax.text(row["logfoldchange"], row["neg_log10_pval"], f" {row['gene']}", fontsize=8)

    ax.set_title(title, fontweight="bold")
    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10 Adjusted P-value")

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight", dpi=300)
    return fig
