"""
Quality Control and Doublet Detection Helper Module.
Provides QC metric calculation, filtering functions, and Scrublet doublet detection.
Complies with scRNAseq_Workbench_Requirements.md (FR-6, NFR-4, NFR-5).
"""

from typing import Optional, Tuple, Union

import anndata as ad
import scanpy as sc
from rich.console import Console

console = Console()


def calculate_qc_metrics(
    adata: ad.AnnData,
    mito_prefix: Union[str, Tuple[str, ...]] = ("MT-", "mt-"),
    ribo_prefix: Union[str, Tuple[str, ...]] = ("RPS", "RPL", "rps", "rpl"),
) -> ad.AnnData:
    """Calculate per-cell and per-gene quality control metrics using Scanpy.

    Adds to adata.obs:
      - n_genes_by_counts : Detected genes count (> 0 counts)
      - total_counts      : Total sequenced UMI / read depth
      - pct_counts_mito   : Percentage of mitochondrial gene counts
      - pct_counts_ribo   : Percentage of ribosomal protein gene counts
    """
    if isinstance(mito_prefix, str):
        mito_prefix = (mito_prefix,)
    if isinstance(ribo_prefix, str):
        ribo_prefix = (ribo_prefix,)

    # Identify mito / ribo genes
    gene_names = adata.var.get("gene_name", adata.var_names)
    adata.var["mt"] = [any(str(g).startswith(p) for p in mito_prefix) for g in gene_names]
    adata.var["ribo"] = [any(str(g).startswith(p) for p in ribo_prefix) for g in gene_names]

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        layer="counts" if "counts" in adata.layers else None,
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    # Standardize column naming
    if "pct_counts_mt" in adata.obs:
        adata.obs["pct_counts_mito"] = adata.obs["pct_counts_mt"]
    if "pct_counts_ribo" in adata.obs:
        adata.obs["pct_counts_ribo"] = adata.obs["pct_counts_ribo"]

    console.print(
        f"[bold blue]QC Calculated:[/bold blue] {adata.n_obs} cells, "
        f"{adata.var['mt'].sum()} mitochondrial genes, "
        f"{adata.var['ribo'].sum()} ribosomal genes."
    )
    return adata


def filter_cells(
    adata: ad.AnnData,
    min_genes: int = 200,
    max_genes: Optional[int] = None,
    min_counts: int = 500,
    max_counts: Optional[int] = None,
    max_pct_mito: float = 20.0,
    max_doublet_score: Optional[float] = None,
) -> ad.AnnData:
    """Filter low-quality cells based on QC thresholds and return a filtered AnnData copy."""
    n_initial = adata.n_obs
    mask = (adata.obs["n_genes_by_counts"] >= min_genes) & (adata.obs["total_counts"] >= min_counts)

    if max_genes is not None:
        mask &= adata.obs["n_genes_by_counts"] <= max_genes

    if max_counts is not None:
        mask &= adata.obs["total_counts"] <= max_counts

    if "pct_counts_mito" in adata.obs:
        mask &= adata.obs["pct_counts_mito"] <= max_pct_mito

    if max_doublet_score is not None and "doublet_score" in adata.obs:
        mask &= adata.obs["doublet_score"] <= max_doublet_score

    filtered_adata = adata[mask].copy()
    n_filtered = filtered_adata.n_obs
    console.print(
        f"[bold green]Cell Filtering:[/bold green] Retained {n_filtered}/{n_initial} cells "
        f"({n_initial - n_filtered} cells removed)."
    )
    return filtered_adata


def filter_genes(adata: ad.AnnData, min_cells: int = 3) -> ad.AnnData:
    """Filter out genes detected in fewer than min_cells."""
    n_initial = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=min_cells)
    console.print(
        f"[bold green]Gene Filtering:[/bold green] Retained {adata.n_vars}/{n_initial} genes "
        f"(min_cells={min_cells})."
    )
    return adata


def detect_doublets(
    adata: ad.AnnData,
    expected_doublet_rate: float = 0.06,
    random_state: int = 42,
    n_prin_comps: Optional[int] = None,
) -> ad.AnnData:
    """Run Scrublet doublet detection directly on raw count matrix.

    Adds to adata.obs:
      - doublet_score    : Continuous doublet score [0.0, 1.0]
      - predicted_doublet: Boolean doublet classification
    """
    import scrublet as scr

    counts_matrix = adata.layers["counts"] if "counts" in adata.layers else adata.X
    console.print(
        f"[bold cyan]Running Scrublet doublet prediction (expected rate: {expected_doublet_rate})...[/bold cyan]"
    )

    scrub = scr.Scrublet(
        counts_matrix,
        expected_doublet_rate=expected_doublet_rate,
        random_state=random_state,
    )

    n_samples, n_vars = counts_matrix.shape
    if n_prin_comps is None:
        # Default is 30 for real single-cell datasets, scaled down for smaller test matrices
        n_pcs = min(30, max(2, (min(n_samples, n_vars) // 6)))
    else:
        n_pcs = n_prin_comps

    doublet_scores, predicted_doublets = scrub.scrub_doublets(n_prin_comps=n_pcs, verbose=False)

    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = (
        predicted_doublets if predicted_doublets is not None else (doublet_scores > 0.3)
    )

    n_doublets = int(adata.obs["predicted_doublet"].sum())
    console.print(
        f"[bold green]✔ Scrublet finished:[/bold green] {n_doublets}/{adata.n_obs} predicted doublets."
    )
    return adata
