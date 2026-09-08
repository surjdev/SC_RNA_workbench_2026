"""
Dimensionality Reduction and Batch Integration Module.
Implements PCA, Harmony batch alignment, k-NN graph construction, UMAP, and t-SNE.
"""

from typing import Optional

import anndata as ad
import scanpy as sc
from rich.console import Console

console = Console()


def run_pca(
    adata: ad.AnnData, n_comps: int = 50, svd_solver: str = "arpack", zero_center: bool = True
) -> ad.AnnData:
    """Compute Principal Component Analysis on highly variable genes."""
    n_features = (
        int(adata.var["highly_variable"].sum()) if "highly_variable" in adata.var else adata.n_vars
    )
    n_components = min(n_comps, adata.n_obs - 1, n_features - 1)
    if n_components < 1:
        raise ValueError("PCA needs at least two cells and two selected genes")
    console.print(f"[bold cyan]Computing PCA ({n_components} components)...[/bold cyan]")

    mask_var = (
        "highly_variable"
        if ("highly_variable" in adata.var and adata.var["highly_variable"].sum() > 0)
        else None
    )
    sc.pp.pca(
        adata,
        n_comps=n_components,
        mask_var=mask_var,
        svd_solver=svd_solver,
        zero_center=zero_center,
    )
    console.print("[bold green]✔ PCA calculation complete (stored in .obsm['X_pca']).[/bold green]")
    return adata


def run_harmony(
    adata: ad.AnnData, batch_key: str = "batch", max_iter_harmony: int = 20
) -> ad.AnnData:
    """
    Run Harmony multi-dataset / multi-batch integration algorithm.
    Saves integrated latent space in .obsm['X_pca_harmony'].
    """
    import harmonypy as hm

    if batch_key not in adata.obs:
        raise KeyError(f"Batch key '{batch_key}' not found in adata.obs!")

    console.print(f"[bold cyan]Running Harmony batch integration on '{batch_key}'...[/bold cyan]")

    if "X_pca" not in adata.obsm:
        adata = run_pca(adata)

    pca_mat = adata.obsm["X_pca"]

    ho = hm.run_harmony(
        pca_mat, adata.obs, batch_key, max_iter_harmony=max_iter_harmony, verbose=False
    )

    adata.obsm["X_pca_harmony"] = ho.Z_corr.T
    console.print(
        "[bold green]✔ Harmony integration complete (stored in .obsm['X_pca_harmony']).[/bold green]"
    )
    return adata


def compute_neighbors(
    adata: ad.AnnData, n_neighbors: int = 15, n_pcs: int = 30, use_rep: Optional[str] = None
) -> ad.AnnData:
    """
    Construct cell-cell k-nearest neighbor (k-NN) neighborhood graph.
    """
    if use_rep is None:
        if "X_pca_harmony" in adata.obsm:
            use_rep = "X_pca_harmony"
        elif "X_pca" in adata.obsm:
            use_rep = "X_pca"
        else:
            adata = run_pca(adata)
            use_rep = "X_pca"

    n_components = min(n_pcs, adata.obsm[use_rep].shape[1])
    console.print(
        f"[bold cyan]Building k-NN graph ({n_neighbors} neighbors, {n_components} components from {use_rep})...[/bold cyan]"
    )

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_components, use_rep=use_rep)
    console.print("[bold green]✔ k-NN neighborhood graph constructed.[/bold green]")
    return adata


def run_umap(
    adata: ad.AnnData, min_dist: float = 0.3, spread: float = 1.0, random_state: int = 42
) -> ad.AnnData:
    """Compute 2D Uniform Manifold Approximation and Projection (UMAP) embedding."""
    if "neighbors" not in adata.uns:
        adata = compute_neighbors(adata)

    console.print("[bold cyan]Computing 2D UMAP non-linear manifold projection...[/bold cyan]")
    sc.tl.umap(adata, min_dist=min_dist, spread=spread, random_state=random_state)
    console.print(
        "[bold green]✔ UMAP coordinates computed (stored in .obsm['X_umap']).[/bold green]"
    )
    return adata


def run_tsne(adata: ad.AnnData, perplexity: float = 30.0, random_state: int = 42) -> ad.AnnData:
    """Compute 2D t-Distributed Stochastic Neighbor Embedding (t-SNE)."""
    console.print("[bold cyan]Computing 2D t-SNE projection...[/bold cyan]")
    sc.tl.tsne(adata, perplexity=perplexity, random_state=random_state)
    console.print(
        "[bold green]✔ t-SNE coordinates computed (stored in .obsm['X_tsne']).[/bold green]"
    )
    return adata
