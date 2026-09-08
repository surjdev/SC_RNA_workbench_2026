"""
Cell Clustering and Community Detection Module.
Provides Leiden and Louvain graph clustering algorithms and sub-clustering utilities.
"""

from typing import Optional, Union

import anndata as ad
import pandas as pd
import scanpy as sc
from rich.console import Console

console = Console()


def cluster_leiden(
    adata: ad.AnnData, resolution: float = 0.6, key_added: str = "leiden", random_state: int = 42
) -> ad.AnnData:
    """
    Cluster cells using the Leiden community detection algorithm.
    """
    if "neighbors" not in adata.uns:
        from .reduction import compute_neighbors

        adata = compute_neighbors(adata)

    console.print(
        f"[bold cyan]Clustering cells with Leiden (resolution={resolution})...[/bold cyan]"
    )
    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added=key_added,
        random_state=random_state,
        flavor="igraph",
        n_iterations=2,
    )

    n_clusters = len(adata.obs[key_added].unique())
    console.print(
        f"[bold green]✔ Leiden clustering complete:[/bold green] Discovered [magenta]{n_clusters}[/magenta] clusters."
    )
    return adata


def cluster_louvain(
    adata: ad.AnnData, resolution: float = 0.6, key_added: str = "louvain", random_state: int = 42
) -> ad.AnnData:
    """
    Cluster cells using the Louvain community detection algorithm.
    """
    if "neighbors" not in adata.uns:
        from .reduction import compute_neighbors

        adata = compute_neighbors(adata)

    console.print(
        f"[bold cyan]Clustering cells with Louvain (resolution={resolution})...[/bold cyan]"
    )
    sc.tl.louvain(adata, resolution=resolution, key_added=key_added, random_state=random_state)

    n_clusters = len(adata.obs[key_added].unique())
    console.print(
        f"[bold green]✔ Louvain clustering complete:[/bold green] Discovered [magenta]{n_clusters}[/magenta] clusters."
    )
    return adata


def subcluster(
    adata: ad.AnnData,
    parent_cluster: Union[str, int],
    cluster_key: str = "leiden",
    resolution: float = 0.4,
    new_key: Optional[str] = None,
) -> ad.AnnData:
    """
    Perform high-resolution sub-clustering on a specific cell subpopulation.
    """
    if new_key is None:
        new_key = f"{cluster_key}_sub_{parent_cluster}"

    console.print(
        f"[bold cyan]Sub-clustering cluster '{parent_cluster}' (resolution={resolution})...[/bold cyan]"
    )

    mask = adata.obs[cluster_key].astype(str) == str(parent_cluster)
    sub_adata = adata[mask].copy()

    sc.pp.neighbors(sub_adata, n_neighbors=10, n_pcs=15)
    sc.tl.leiden(sub_adata, resolution=resolution, key_added=new_key)

    # Map back to parent AnnData
    adata.obs[new_key] = adata.obs[cluster_key].astype(str)
    adata.obs.loc[mask, new_key] = [f"{parent_cluster}.{c}" for c in sub_adata.obs[new_key]]

    console.print(
        f"[bold green]✔ Sub-clustering complete (saved in .obs['{new_key}']).[/bold green]"
    )
    return adata


def cluster_sklearn(
    adata: ad.AnnData, estimator, use_rep: str = "X_pca", key_added: str = "sklearn_cluster"
) -> ad.AnnData:
    """
    Cluster cells using any Scikit-Learn clustering algorithm or estimator.

    Parameters
    ----------
    adata : ad.AnnData
        Annotated dataset.
    estimator : sklearn.base.ClusterMixin or object with fit_predict / fit
        Scikit-learn clustering instance (e.g. KMeans, AgglomerativeClustering, DBSCAN).
    use_rep : str, default 'X_pca'
        Representation key in adata.obsm to cluster on (e.g. 'X_pca', 'X_umap', 'X_pca_harmony').
    key_added : str, default 'sklearn_cluster'
        Column name in adata.obs to store assigned cluster labels.
    """
    import numpy as np
    import scipy.sparse as sp

    if use_rep in adata.obsm:
        X = adata.obsm[use_rep]
    elif use_rep == "X":
        X = adata.X.toarray() if sp.issparse(adata.X) else adata.X
    else:
        from .reduction import run_pca

        adata = run_pca(adata)
        X = adata.obsm["X_pca"]

    est_name = estimator.__class__.__name__
    console.print(
        f"[bold cyan]Clustering cells using Scikit-Learn [yellow]{est_name}[/yellow] on [magenta]{use_rep}[/magenta]...[/bold cyan]"
    )

    if hasattr(estimator, "fit_predict"):
        labels = estimator.fit_predict(X)
    else:
        estimator.fit(X)
        labels = getattr(estimator, "labels_", None)
        if labels is None:
            raise ValueError(f"Estimator {est_name} has no labels_ attribute after fit()!")

    adata.obs[key_added] = pd.Categorical([str(lbl) for lbl in labels])
    n_clusters = len(np.unique(labels))
    console.print(
        f"[bold green]✔ Scikit-Learn ({est_name}) clustering complete:[/bold green] Discovered [magenta]{n_clusters}[/magenta] clusters (stored in .obs['{key_added}'])."
    )
    return adata


def cluster_kmeans(
    adata: ad.AnnData,
    n_clusters: int = 5,
    use_rep: str = "X_pca",
    key_added: str = "kmeans",
    random_state: int = 42,
) -> ad.AnnData:
    """Cluster cells using Scikit-Learn KMeans."""
    from sklearn.cluster import KMeans

    km = KMeans(n_clusters=n_clusters, random_state=random_state, n_init="auto")
    return cluster_sklearn(adata, estimator=km, use_rep=use_rep, key_added=key_added)


def cluster_hierarchical(
    adata: ad.AnnData, n_clusters: int = 5, use_rep: str = "X_pca", key_added: str = "hierarchical"
) -> ad.AnnData:
    """Cluster cells using Scikit-Learn AgglomerativeClustering."""
    from sklearn.cluster import AgglomerativeClustering

    agg = AgglomerativeClustering(n_clusters=n_clusters)
    return cluster_sklearn(adata, estimator=agg, use_rep=use_rep, key_added=key_added)
