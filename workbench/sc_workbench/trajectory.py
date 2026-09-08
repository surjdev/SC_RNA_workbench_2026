"""
Lineage Trajectory and Pseudotime Inference Module.
Implements Partition-based Graph Abstraction (PAGA) and Diffusion Pseudotime (DPT).
"""

from typing import Optional, Union
import anndata as ad
import scanpy as sc
from rich.console import Console

console = Console()

def run_paga(
    adata: ad.AnnData,
    groups: str = "leiden",
    use_rna_velocity: bool = False
) -> ad.AnnData:
    """
    Compute PAGA (Partition-based Graph Abstraction) to infer cellular lineage topology.
    """
    if "neighbors" not in adata.uns:
        from .reduction import compute_neighbors
        adata = compute_neighbors(adata)
        
    console.print(f"[bold cyan]Running PAGA graph abstraction on '{groups}'...[/bold cyan]")
    sc.tl.paga(
        adata,
        groups=groups,
        use_rna_velocity=use_rna_velocity
    )
    console.print("[bold green]✔ PAGA connectivity matrix computed (stored in .uns['paga']).[/bold green]")
    return adata


def run_dpt(
    adata: ad.AnnData,
    root_cell: Optional[str] = None,
    root_cluster: Optional[Union[str, int]] = None,
    cluster_key: str = "leiden",
    n_dcs: int = 10
) -> ad.AnnData:
    """
    Calculate Diffusion Pseudotime (DPT) from an inferred or user-defined progenitor root cell.
    Adds `dpt_pseudotime` to .obs.
    """
    console.print("[bold cyan]Computing Diffusion Maps and Diffusion Pseudotime (DPT)...[/bold cyan]")
    
    if root_cell is None and root_cluster is None:
        raise ValueError("Provide root_cell or root_cluster based on biological evidence")
    if root_cell is not None and root_cluster is not None:
        raise ValueError("Choose only one of root_cell or root_cluster")
    # 1. Compute Diffusion Map
    sc.tl.diffmap(adata, n_comps=n_dcs)
    
    # 2. Select root cell
    if root_cell is not None:
        if root_cell not in adata.obs_names:
            raise KeyError(f"Root cell '{root_cell}' not found in adata.obs_names!")
        adata.uns["iroot"] = adata.obs_names.get_loc(root_cell)
    elif root_cluster is not None:
        cluster_cells = adata.obs[adata.obs[cluster_key].astype(str) == str(root_cluster)].index
        if len(cluster_cells) == 0:
            raise ValueError(f"No cells found in root cluster '{root_cluster}'!")
        adata.uns["iroot"] = adata.obs_names.get_loc(cluster_cells[0])
    else:
        # Automatically select extreme cell along DC1
        adata.uns["iroot"] = int(adata.obsm["X_diffmap"][:, 1].argmin())
        
    # 3. Compute DPT
    sc.tl.dpt(adata, n_dcs=n_dcs)
    console.print("[bold green]✔ Diffusion pseudotime computed (stored in .obs['dpt_pseudotime']).[/bold green]")
    return adata
