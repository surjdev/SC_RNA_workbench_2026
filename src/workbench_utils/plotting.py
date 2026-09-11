"""Publication-Quality Plotting and Visualization Module for SMART-seq2.

Provides Nature/Cell journal color schemes, QC violin plots with ERCC spike-ins,
plate well layout visualizations, and publication-ready Volcano plots.
Complies with scRNAseq_Workbench_Requirements.md (FR-6, NFR-4, NFR-5).
"""

from typing import Dict, Optional, Tuple

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

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
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
]


def set_publication_style(font_scale: float = 1.0, style: str = "whitegrid") -> None:
    """Set global matplotlib and seaborn publication styling."""
    sns.set_theme(style=style, font_scale=font_scale)
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Helvetica", "Arial"]
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["figure.dpi"] = 150
    plt.rcParams["axes.edgecolor"] = "#333333"
    plt.rcParams["axes.linewidth"] = 0.8
    plt.rcParams["xtick.major.width"] = 0.8
    plt.rcParams["ytick.major.width"] = 0.8


def plot_qc_violins(
    adata: ad.AnnData,
    keys: Tuple[str, ...] = (
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mito",
        "pct_counts_ercc",
    ),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Plot publication-quality violin distributions for SMART-seq2 QC metrics."""
    set_publication_style()

    available_keys = [k for k in keys if k in adata.obs]
    if not available_keys:
        raise KeyError(f"None of {keys} found in adata.obs!")

    fig, axes = plt.subplots(1, len(available_keys), figsize=(4 * len(available_keys), 4.2))
    if len(available_keys) == 1:
        axes = [axes]

    titles: Dict[str, str] = {
        "n_genes_by_counts": "Detected Genes per Cell",
        "total_counts": "Total Read Counts",
        "pct_counts_mito": "Mitochondrial Reads (%)",
        "pct_counts_ribo": "Ribosomal Reads (%)",
        "pct_counts_ercc": "ERCC Spike-ins (%)",
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


def plot_plate_layout(
    adata: ad.AnnData,
    plate_id: Optional[str] = None,
    color_key: str = "total_counts",
    plate_key: str = "plate",
    title: Optional[str] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Visualize 96-well plate layout (Rows A-H, Cols 1-12) colored by a QC metric."""
    set_publication_style()

    sub_adata = adata
    if plate_id and plate_key in adata.obs:
        sub_adata = adata[adata.obs[plate_key] == plate_id]

    rows = ["A", "B", "C", "D", "E", "F", "G", "H"]
    cols = list(range(1, 13))

    grid = pd.DataFrame(np.nan, index=rows, columns=cols)

    for _, cell in sub_adata.obs.iterrows():
        r = cell.get("well_row")
        c = cell.get("well_col")
        val = cell.get(color_key)
        if pd.notna(r) and pd.notna(c) and r in rows and int(c) in cols:
            grid.loc[r, int(c)] = val

    fig, ax = plt.subplots(figsize=(9, 5))
    has_valid_values = np.any(np.isfinite(grid.values))
    max_val = np.nanmax(grid.values) if has_valid_values else 0.0

    sns.heatmap(
        grid,
        annot=True,
        fmt=".0f" if max_val > 50 else ".1f",
        cmap="YlGnBu",
        cbar_kws={"label": color_key},
        ax=ax,
        linewidths=0.5,
        linecolor="#cccccc",
    )
    chart_title = (
        title if title is not None else f"Plate Layout: {plate_id or 'Plate 1'} ({color_key})"
    )
    ax.set_title(chart_title, fontsize=13, fontweight="bold")
    ax.set_xlabel("Plate Column (1-12)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Plate Row (A-H)", fontsize=11, fontweight="bold")

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
    """Generate publication-ready Volcano Plot for PyDESeq2 / DE results."""
    set_publication_style()
    fig, ax = plt.subplots(figsize=(7, 6))

    df = de_df.copy()
    if gene_col not in df.columns:
        df[gene_col] = df.index

    # Drop invalid values
    df = df.dropna(subset=[pval_col, fc_col])
    # Avoid log(0)
    min_pval = df[df[pval_col] > 0][pval_col].min() if (df[pval_col] > 0).any() else 1e-300
    df["-log10(padj)"] = -np.log10(df[pval_col].replace(0, min_pval * 0.1))

    # Determine significance categories
    df["status"] = "Not Significant"
    df.loc[(df[pval_col] < pval_threshold) & (df[fc_col] >= fc_threshold), "status"] = "Upregulated"
    df.loc[(df[pval_col] < pval_threshold) & (df[fc_col] <= -fc_threshold), "status"] = (
        "Downregulated"
    )

    color_map = {
        "Not Significant": "#999999",
        "Upregulated": "#E64B35",
        "Downregulated": "#4DBBD5",
    }

    for status, grp in df.groupby("status"):
        ax.scatter(
            grp[fc_col],
            grp["-log10(padj)"],
            c=color_map.get(status, "#999999"),
            label=f"{status} ({len(grp)})",
            alpha=0.75,
            edgecolor="none",
            s=25,
        )

    # Threshold lines
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
