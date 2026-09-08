"""
Preprocessing and Normalization Module for Single-Cell Analysis.
Handles library depth scaling, log1p transformation, Highly Variable Gene (HVG) selection,
and Z-score feature scaling.
"""

from typing import Optional

import anndata as ad
import scanpy as sc
from rich.console import Console

console = Console()


def normalize_and_log(
    adata: ad.AnnData, target_sum: float = 1e4, save_raw: bool = True
) -> ad.AnnData:
    """
    Perform depth normalization and natural log transformation.
    Saves the pre-scaling normalized data in .raw or .layers['normalized'].
    """
    console.print(
        f"[bold cyan]Normalizing counts to {target_sum:,.0f} and computing log1p...[/bold cyan]"
    )

    # Store raw counts in layer if not present
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    adata.uns.pop("log1p", None)
    sc.pp.log1p(adata)
    adata.layers["normalized"] = adata.X.copy()

    if save_raw:
        adata.raw = adata

    console.print("[bold green]✔ Normalization & log1p complete.[/bold green]")
    return adata


def select_hvg(
    adata: ad.AnnData,
    n_top_genes: int = 2000,
    flavor: str = "seurat",
    batch_key: Optional[str] = None,
) -> ad.AnnData:
    """
    Identify Highly Variable Genes (HVGs) for dimensionality reduction.
    """
    console.print(
        f"[bold cyan]Selecting top {n_top_genes} Highly Variable Genes (flavor={flavor})...[/bold cyan]"
    )

    # Check if number of genes is smaller than n_top_genes
    n_genes = min(n_top_genes, adata.n_vars)

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_genes,
        flavor=flavor,
        layer="counts" if flavor in ("seurat_v3", "seurat_v3_paper") else None,
        batch_key=batch_key,
        subset=False,
    )

    n_hvg = int(adata.var["highly_variable"].sum())
    console.print(
        f"[bold green]✔ HVG Selection:[/bold green] Identified [yellow]{n_hvg}[/yellow] highly variable genes."
    )
    return adata


def scale_features(
    adata: ad.AnnData, max_value: float = 10.0, zero_center: bool = True
) -> ad.AnnData:
    """
    Standardize expression of genes to unit variance and zero mean.
    Saves scaled matrix in .layers['scaled'] and updates .X.
    """
    console.print("[bold cyan]Standardizing gene features (Z-score scaling)...[/bold cyan]")
    sc.pp.scale(adata, max_value=max_value, zero_center=zero_center)
    adata.layers["scaled"] = adata.X.copy()
    console.print("[bold green]✔ Feature scaling complete.[/bold green]")
    return adata
