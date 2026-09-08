"""
Marker Gene Discovery and Differential Expression Module.
Ranks cluster-specific biomarkers using Wilcoxon rank-sum tests, Welch t-tests,
and extracts tidy DataFrames for downstream reporting.
"""

from typing import Optional

import anndata as ad
import pandas as pd
import scanpy as sc
from rich.console import Console

console = Console()


def find_markers(
    adata: ad.AnnData,
    groupby: str = "leiden",
    method: str = "wilcoxon",
    n_genes: int = 50,
    pts: bool = True,
    key_added: str = "rank_genes_groups",
) -> ad.AnnData:
    """
    Identify cluster biomarker genes by computing one-vs-rest differential expression.
    """
    console.print(
        f"[bold cyan]Finding cluster marker genes using {method} test on '{groupby}'...[/bold cyan]"
    )

    # Use normalized/log1p layer or .raw if available
    use_raw = adata.raw is not None

    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        n_genes=n_genes,
        pts=pts,
        key_added=key_added,
        layer="normalized" if not use_raw and "normalized" in adata.layers else None,
        use_raw=use_raw,
    )

    console.print("[bold green]✔ Marker gene discovery complete.[/bold green]")
    return adata


def get_markers_df(
    adata: ad.AnnData,
    key: str = "rank_genes_groups",
    pval_cutoff: Optional[float] = None,
    logfc_cutoff: Optional[float] = None,
) -> pd.DataFrame:
    """
    Extract differential expression results from AnnData into a tidy Pandas DataFrame.
    """
    if key not in adata.uns:
        raise KeyError(f"Key '{key}' not found in adata.uns! Run find_markers() first.")

    result = adata.uns[key]
    groups = result["names"].dtype.names

    records = []
    for group in groups:
        for idx in range(len(result["names"][group])):
            gene = result["names"][group][idx]
            logfc = result["logfoldchanges"][group][idx]
            pval = result["pvals"][group][idx]
            pval_adj = result["pvals_adj"][group][idx]
            score = result["scores"][group][idx]

            # Percent expressing if available
            pct_group = result["pts"][group].loc[gene] if "pts" in result else None
            pct_rest = result["pts_rest"][group].loc[gene] if "pts_rest" in result else None

            keep = True
            if pval_cutoff is not None and pval_adj > pval_cutoff:
                keep = False
            if logfc_cutoff is not None and abs(logfc) < logfc_cutoff:
                keep = False

            if keep:
                records.append(
                    {
                        "cluster": str(group),
                        "gene": str(gene),
                        "score": float(score),
                        "logfoldchange": float(logfc),
                        "pvals": float(pval),
                        "pvals_adj": float(pval_adj),
                        "pct_nz_group": float(pct_group) if pct_group is not None else None,
                        "pct_nz_reference": float(pct_rest) if pct_rest is not None else None,
                    }
                )

    df = pd.DataFrame(records)
    if df.empty:
        df = pd.DataFrame(
            columns=[
                "cluster",
                "gene",
                "score",
                "logfoldchange",
                "pvals",
                "pvals_adj",
                "pct_nz_group",
                "pct_nz_reference",
            ]
        )
    return df


def find_pairwise_markers(
    adata: ad.AnnData, group1: str, group2: str, groupby: str = "leiden", method: str = "wilcoxon"
) -> pd.DataFrame:
    """
    Compute differential expression directly comparing group1 vs group2.
    """
    console.print(f"[bold cyan]Comparing '{group1}' vs '{group2}' on '{groupby}'...[/bold cyan]")
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        groups=[group1],
        reference=group2,
        method=method,
        key_added="pairwise_de",
    )
    df = get_markers_df(adata, key="pairwise_de", pval_cutoff=1.0, logfc_cutoff=0.0)
    return df
