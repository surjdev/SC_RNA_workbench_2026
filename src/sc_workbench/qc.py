"""
Quality Control and Doublet Detection Module for Single-Cell Analysis.
Calculates library depth, detected gene count, mitochondrial/ribosomal fractions,
and performs automated doublet prediction using Scrublet.
"""

from typing import Optional, Tuple, Union

import anndata as ad
import numpy as np
import scanpy as sc
from rich.console import Console

console = Console()


def calculate_qc_metrics(
    adata: ad.AnnData,
    mito_prefix: Union[str, Tuple[str, ...]] = ("MT-", "mt-"),
    ribo_prefix: Union[str, Tuple[str, ...]] = ("RPS", "RPL", "rps", "rpl"),
) -> ad.AnnData:
    """
    Calculate per-cell and per-gene quality control metrics.
    Adds to .obs:
      - n_genes_by_counts : Number of genes with > 0 counts
      - total_counts      : Total sequenced read / UMI counts
      - pct_counts_mito   : Percentage of reads mapped to mitochondrial genes
      - pct_counts_ribo   : Percentage of reads mapped to ribosomal genes
    """
    if isinstance(mito_prefix, str):
        mito_prefix = (mito_prefix,)
    if isinstance(ribo_prefix, str):
        ribo_prefix = (ribo_prefix,)

    # Tag mitochondrial genes
    adata.var["mt"] = [
        any(str(g).startswith(p) for p in mito_prefix)
        for g in adata.var.get("gene_name", adata.var_names)
    ]
    # Tag ribosomal protein genes
    adata.var["ribo"] = [
        any(str(g).startswith(p) for p in ribo_prefix)
        for g in adata.var.get("gene_name", adata.var_names)
    ]

    # Calculate Scanpy standard QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        layer="counts" if "counts" in adata.layers else None,
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    # Rename for cleaner standard naming
    if "pct_counts_mt" in adata.obs:
        adata.obs["pct_counts_mito"] = adata.obs["pct_counts_mt"]
    if "pct_counts_ribo" in adata.obs:
        adata.obs["pct_counts_ribo"] = adata.obs["pct_counts_ribo"]

    console.print(
        f"[bold blue]QC Calculated:[/bold blue] "
        f"Median counts: [green]{int(np.median(adata.obs['total_counts'])):,}[/green] | "
        f"Median genes: [green]{int(np.median(adata.obs['n_genes_by_counts'])):,}[/green] | "
        f"Median Mito %: [green]{np.median(adata.obs['pct_counts_mito']):.2f}%[/green]"
    )
    return adata


def detect_doublets(
    adata: ad.AnnData, expected_doublet_rate: float = 0.06, random_state: int = 42
) -> ad.AnnData:
    """
    Predict cellular doublets using the Scrublet algorithm.
    Adds to .obs:
      - doublet_score : Continuous doublet likelihood score [0.0 - 1.0]
      - is_doublet    : Boolean classification flag
    """
    import scrublet as scr

    console.print("[bold cyan]Running Scrublet doublet detection...[/bold cyan]")
    try:
        # Use raw counts layer if available
        counts_matrix = adata.layers["counts"] if "counts" in adata.layers else adata.X

        scrub = scr.Scrublet(
            counts_matrix, expected_doublet_rate=expected_doublet_rate, random_state=random_state
        )
        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30, verbose=False
        )

        if predicted_doublets is None:
            # Threshold fallback if automated bimodality check was flat
            raise RuntimeError("Scrublet could not infer a threshold; inspect scores explicitly")

        adata.obs["doublet_score"] = doublet_scores
        adata.obs["is_doublet"] = predicted_doublets.astype(bool)

        n_doublets = int(adata.obs["is_doublet"].sum())
        console.print(
            f"[bold green]✔ Doublets detected:[/bold green] {n_doublets} / {adata.n_obs} cells ({n_doublets / adata.n_obs * 100:.1f}%)"
        )
    except Exception as e:
        raise RuntimeError("Scrublet failed; no doublet calls were assigned") from e

    return adata


def filter_cells(
    adata: ad.AnnData,
    min_genes: int = 200,
    max_genes: Optional[int] = None,
    min_counts: int = 500,
    max_counts: Optional[int] = None,
    max_pct_mito: float = 20.0,
    filter_doublets: bool = False,
) -> ad.AnnData:
    """
    Filter low-quality cells based on gene counts, total read depth, mitochondrial percentage,
    and doublet calls.
    """
    initial_cells = adata.n_obs

    if "total_counts" not in adata.obs or "n_genes_by_counts" not in adata.obs:
        adata = calculate_qc_metrics(adata)

    # Build filtering mask
    mask = (
        (adata.obs["n_genes_by_counts"] >= min_genes)
        & (adata.obs["total_counts"] >= min_counts)
        & (adata.obs["pct_counts_mito"] <= max_pct_mito)
    )

    if max_genes:
        mask = mask & (adata.obs["n_genes_by_counts"] <= max_genes)
    if max_counts:
        mask = mask & (adata.obs["total_counts"] <= max_counts)
    if filter_doublets and "is_doublet" not in adata.obs:
        raise ValueError("Run doublet detection before filtering doublets")
    if filter_doublets and "is_doublet" in adata.obs:
        mask = mask & (~adata.obs["is_doublet"])

    adata_filtered = adata[mask].copy()
    retained_cells = adata_filtered.n_obs
    dropped_cells = initial_cells - retained_cells

    console.print(
        f"[bold green]✔ Filtered Cells:[/bold green] Retained [cyan]{retained_cells}[/cyan] cells "
        f"(dropped [yellow]{dropped_cells}[/yellow] / {initial_cells}, "
        f"[bold]{retained_cells / initial_cells * 100:.1f}% retained[/bold])"
    )
    return adata_filtered


def filter_genes(adata: ad.AnnData, min_cells: int = 3) -> ad.AnnData:
    """Filter out genes that are detected in fewer than `min_cells` cells."""
    initial_genes = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=min_cells)
    console.print(
        f"[bold green]✔ Filtered Genes:[/bold green] Retained {adata.n_vars} / {initial_genes} genes (min_cells >= {min_cells})"
    )
    return adata
