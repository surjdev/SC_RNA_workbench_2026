"""
Cell Type Annotation and Signature Scoring Module.
Enables signature activity scoring across known cell markers and mapping cluster identities
to curated biological cell type nomenclatures.
"""

from typing import Dict, List, Union

import anndata as ad
import scanpy as sc
from rich.console import Console

console = Console()


def score_gene_signatures(
    adata: ad.AnnData, signatures: Dict[str, List[str]], ctrl_size: int = 50, prefix: str = "sig_"
) -> ad.AnnData:
    """
    Score cell-specific activity for sets of marker genes.
    Adds scores to .obs under `<prefix><signature_name>`.
    """
    console.print(
        f"[bold cyan]Scoring {len(signatures)} gene signatures across cells...[/bold cyan]"
    )

    for sig_name, genes in signatures.items():
        # Keep genes present in adata
        valid_genes = [g for g in genes if g in adata.var_names]
        if len(valid_genes) == 0:
            console.print(
                f"[yellow]⚠ No valid genes found in dataset for signature '{sig_name}'![/yellow]"
            )
            continue

        score_key = f"{prefix}{sig_name}"
        sc.tl.score_genes(
            adata,
            gene_list=valid_genes,
            score_name=score_key,
            ctrl_size=min(ctrl_size, len(adata.var_names) - len(valid_genes)),
            use_raw=(adata.raw is not None),
        )
        console.print(
            f"[bold green]✔ Scored '{sig_name}':[/bold green] {len(valid_genes)} genes -> .obs['{score_key}']"
        )

    return adata


def assign_cell_types(
    adata: ad.AnnData,
    cluster_to_celltype: Dict[Union[str, int], str],
    cluster_key: str = "leiden",
    new_key: str = "cell_type",
) -> ad.AnnData:
    """
    Map cluster numerical or string identifiers to biological cell types.
    """
    if cluster_key not in adata.obs:
        raise KeyError(f"Cluster key '{cluster_key}' not found in adata.obs!")

    console.print(f"[bold cyan]Mapping clusters in '{cluster_key}' to '{new_key}'...[/bold cyan]")

    # Create mapping string
    str_map = {str(k): v for k, v in cluster_to_celltype.items()}

    adata.obs[new_key] = (
        adata.obs[cluster_key]
        .astype(str)
        .map(str_map)
        .fillna(adata.obs[cluster_key].astype(str))
        .astype("category")
    )

    cell_types = adata.obs[new_key].value_counts()
    console.print("[bold green]✔ Cell types assigned successfully:[/bold green]")
    for ct, count in cell_types.items():
        console.print(
            f"  • [magenta]{ct}[/magenta]: {count} cells ({count / adata.n_obs * 100:.1f}%)"
        )

    return adata
