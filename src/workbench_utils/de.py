"""
Differential Expression Helper Module utilizing PyDESeq2.
Facilitates AnnData conversion to PyDESeq2-ready formats and provides automated DE execution.
Complies with scRNAseq_Workbench_Requirements.md (FR-2, FR-6, NFR-5).
"""

from typing import List, Optional, Tuple, Union

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from rich.console import Console

console = Console()


def prepare_pydeseq2_data(
    adata: ad.AnnData,
    design_factor: str,
    layer: Optional[str] = "counts",
    min_cells_per_gene: int = 3,
    sample_key: Optional[str] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Extract and format raw integer counts and sample metadata for PyDESeq2 DeseqDataSet.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData containing single-cell counts.
    design_factor : str
        Column in adata.obs representing the primary biological comparison variable.
    layer : str, optional
        Layer to extract raw counts from (default is 'counts'). If not found, uses .X.
    min_cells_per_gene : int, default 3
        Pre-filter genes expressed in fewer than this number of cells.
    sample_key : str, optional
        Column in adata.obs to aggregate cells into pseudo-bulk samples.
        If None, each cell is treated as an individual observation.

    Returns
    -------
    counts_df : pd.DataFrame
        Cells/Samples as rows, Genes as columns (integer values).
    clinical_df : pd.DataFrame
        Metadata dataframe indexed matching counts_df.
    """
    if design_factor not in adata.obs:
        raise KeyError(f"Design factor '{design_factor}' not found in adata.obs!")

    # 1. Extract raw counts
    if layer is not None and layer in adata.layers:
        X = adata.layers[layer]
    else:
        X = adata.X

    if sp.issparse(X):
        X = X.toarray()

    # PyDESeq2 requires non-negative integers
    X_int = np.rint(np.asarray(X)).astype(int)

    # 2. Pseudo-bulk aggregation if sample_key provided
    if sample_key is not None:
        if sample_key not in adata.obs:
            raise KeyError(f"Sample key '{sample_key}' not found in adata.obs!")

        console.print(
            f"[bold cyan]Aggregating cells into pseudo-bulk samples by '{sample_key}'...[/bold cyan]"
        )
        sample_ids = adata.obs[sample_key].astype(str)
        unique_samples = np.unique(sample_ids)

        pb_counts = []
        pb_meta = []
        for sid in unique_samples:
            idx = (sample_ids == sid).values
            pb_counts.append(X_int[idx, :].sum(axis=0))
            # Take first observation for metadata factor
            factor_val = adata.obs.loc[idx, design_factor].iloc[0]
            pb_meta.append({sample_key: sid, design_factor: factor_val})

        counts_df = pd.DataFrame(pb_counts, index=unique_samples, columns=adata.var_names)
        clinical_df = pd.DataFrame(pb_meta, index=unique_samples)
    else:
        counts_df = pd.DataFrame(X_int, index=adata.obs_names, columns=adata.var_names)
        clinical_df = adata.obs[[design_factor]].copy()

    # 3. Filter genes with low cell/sample count
    gene_mask = (counts_df > 0).sum(axis=0) >= min_cells_per_gene
    counts_df = counts_df.loc[:, gene_mask]

    console.print(
        f"[bold green]✔ PyDESeq2 data prepared:[/bold green] {counts_df.shape[0]} samples/cells, "
        f"{counts_df.shape[1]} genes (design: '{design_factor}')."
    )
    return counts_df, clinical_df


def run_pydeseq2(
    counts_df: pd.DataFrame,
    clinical_df: pd.DataFrame,
    design_factors: Union[str, List[str]],
    contrast: Optional[Tuple[str, str, str]] = None,
    quiet: bool = True,
) -> pd.DataFrame:
    """Execute PyDESeq2 pipeline and return a tidy DataFrame of differential expression results.

    Parameters
    ----------
    counts_df : pd.DataFrame
        Counts matrix (Samples/Cells x Genes).
    clinical_df : pd.DataFrame
        Metadata dataframe matching counts_df index.
    design_factors : str or list of str
        Factor(s) to include in the design formula.
    contrast : tuple of (str, str, str), optional
        Contrast specification: (factor_name, test_level, reference_level).
    quiet : bool, default True
        Suppress intermediate PyDESeq2 logs.

    Returns
    -------
    pd.DataFrame
        Tidy results containing gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    if isinstance(design_factors, str):
        factors = [design_factors]
    else:
        factors = list(design_factors)

    formula = f"~ {' + '.join(factors)}"
    console.print(f"[bold cyan]Running PyDESeq2 with design formula: {formula}...[/bold cyan]")

    try:
        dds = DeseqDataSet(
            counts=counts_df,
            metadata=clinical_df,
            design=formula,
            quiet=quiet,
        )
    except TypeError:
        dds = DeseqDataSet(
            counts=counts_df,
            metadata=clinical_df,
            design_factors=factors,
            quiet=quiet,
        )
    dds.deseq2()

    stat_res = DeseqStats(dds, contrast=contrast, quiet=quiet)
    stat_res.summary()

    results_df = stat_res.results_df.copy()
    results_df.index.name = "gene"
    results_df = results_df.reset_index()

    # Sort by adjusted p-value
    if "padj" in results_df.columns:
        results_df = results_df.sort_values(by="padj", ascending=True)

    console.print(
        f"[bold green]✔ PyDESeq2 complete:[/bold green] {len(results_df)} genes analyzed."
    )
    return results_df
