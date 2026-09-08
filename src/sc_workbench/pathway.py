"""
Pathway Enrichment and Functional Annotation Module.
Integrates with GSEAPY to perform Over-Representation Analysis (Enrichr)
and Pre-ranked Gene Set Enrichment Analysis (GSEA).
"""

from typing import List, Optional, Union

import pandas as pd
from rich.console import Console

console = Console()


def run_enrichr(
    gene_list: List[str],
    gene_sets: Union[str, List[str]] = "GO_Biological_Process_2023",
    organism: str = "human",
    outdir: Optional[str] = None,
    cutoff: float = 0.05,
    background=None,
) -> pd.DataFrame:
    """
    Perform Over-Representation Analysis (ORA) on a set of marker genes using Enrichr.
    Common libraries:
      - GO_Biological_Process_2023
      - KEGG_2021_Human
      - MSigDB_Hallmark_2020
      - Reactome_2022
    """
    import gseapy as gp

    if isinstance(gene_sets, str):
        gene_sets = [gene_sets]

    console.print(
        f"[bold cyan]Running Enrichr on {len(gene_list)} genes ({gene_sets})...[/bold cyan]"
    )

    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets,
            organism=organism,
            outdir=outdir,
            cutoff=cutoff,
            background=background,
            verbose=False,
        )
        res_df = enr.results
        console.print(
            f"[bold green]✔ Enrichr analysis complete:[/bold green] Found {len(res_df)} terms (unfiltered; cutoff controls plotting)."
        )
        return res_df
    except Exception as e:
        console.print(f"[bold red]✘ Enrichr failed (check internet connection):[/bold red] {e}")
        raise RuntimeError("Enrichment failed; no result was produced") from e


def run_prerank_gsea(
    rnk_series: pd.Series,
    gene_sets: str = "MSigDB_Hallmark_2020",
    min_size: int = 5,
    max_size: int = 500,
    permutation_num: int = 1000,
) -> pd.DataFrame:
    """
    Perform pre-ranked GSEA on gene ranking metrics (e.g. log2 fold change or test statistic).
    """
    import gseapy as gp

    console.print(
        f"[bold cyan]Running Pre-ranked GSEA with {len(rnk_series)} ranked genes...[/bold cyan]"
    )

    try:
        prerank = gp.prerank(
            rnk=rnk_series,
            gene_sets=gene_sets,
            min_size=min_size,
            max_size=max_size,
            permutation_num=permutation_num,
            verbose=False,
            seed=42,
        )
        res_df = prerank.res2d
        console.print(
            f"[bold green]✔ Pre-ranked GSEA complete:[/bold green] Evaluated {len(res_df)} pathways."
        )
        return res_df
    except Exception as e:
        console.print(f"[bold red]✘ GSEA failed:[/bold red] {e}")
        raise RuntimeError("Enrichment failed; no result was produced") from e
