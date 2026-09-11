"""Quality Control and Plate-based Helper Module for SMART-seq2.

Provides QC metric calculation (mitochondrial, ribosomal, ERCC spike-ins),
plate/well coordinate parsing, and cell/gene filtering.
Complies with scRNAseq_Workbench_Requirements.md (FR-6, NFR-4, NFR-5).
"""

import re
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
    ercc_prefix: Union[str, Tuple[str, ...]] = ("ERCC-", "ercc-"),
) -> ad.AnnData:
    """Calculate per-cell and per-gene quality control metrics for SMART-seq2.

    Adds to adata.obs:
      - n_genes_by_counts : Detected genes count (> 0 reads)
      - total_counts      : Total sequenced read depth
      - pct_counts_mito   : Percentage of mitochondrial gene reads
      - pct_counts_ribo   : Percentage of ribosomal protein gene reads
      - pct_counts_ercc   : Percentage of synthetic ERCC spike-in reads
    """
    if isinstance(mito_prefix, str):
        mito_prefix = (mito_prefix,)
    if isinstance(ribo_prefix, str):
        ribo_prefix = (ribo_prefix,)
    if isinstance(ercc_prefix, str):
        ercc_prefix = (ercc_prefix,)

    # Identify gene subsets
    gene_names = adata.var.get("gene_name", adata.var_names)
    adata.var["mt"] = [any(str(g).startswith(p) for p in mito_prefix) for g in gene_names]
    adata.var["ribo"] = [any(str(g).startswith(p) for p in ribo_prefix) for g in gene_names]
    adata.var["ercc"] = [any(str(g).startswith(p) for p in ercc_prefix) for g in gene_names]

    qc_vars = []
    if adata.var["mt"].sum() > 0:
        qc_vars.append("mt")
    if adata.var["ribo"].sum() > 0:
        qc_vars.append("ribo")
    if adata.var["ercc"].sum() > 0:
        qc_vars.append("ercc")

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=qc_vars,
        layer="counts" if "counts" in adata.layers else None,
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    # Standardize column naming
    if "pct_counts_mt" in adata.obs:
        adata.obs["pct_counts_mito"] = adata.obs["pct_counts_mt"]
    else:
        adata.obs["pct_counts_mito"] = 0.0

    if "pct_counts_ribo" in adata.obs:
        adata.obs["pct_counts_ribo"] = adata.obs["pct_counts_ribo"]
    else:
        adata.obs["pct_counts_ribo"] = 0.0

    if "pct_counts_ercc" in adata.obs:
        adata.obs["pct_counts_ercc"] = adata.obs["pct_counts_ercc"]
    else:
        adata.obs["pct_counts_ercc"] = 0.0

    console.print(
        f"[bold blue]SMART-seq2 QC Calculated:[/bold blue] {adata.n_obs} cells, "
        f"{adata.var['mt'].sum()} mito genes, "
        f"{adata.var['ribo'].sum()} ribo genes, "
        f"{adata.var['ercc'].sum()} ERCC spike-ins."
    )
    return adata


def parse_plate_metadata(
    adata: ad.AnnData,
    cell_id_col: Optional[str] = None,
) -> ad.AnnData:
    """Parse plate ID and well coordinates (e.g., 'Plate1_A01', 'P2-H12', 'A01')

    and add 'plate', 'well', 'well_row', 'well_col' columns to adata.obs.
    """
    ids = adata.obs[cell_id_col].astype(str) if cell_id_col else adata.obs_names.astype(str)

    plates = []
    wells = []
    rows = []
    cols = []

    pattern = re.compile(r"^(?:(.*)[_-])?([A-Ha-h])([0-1]?[0-9]|2[0-4])$")

    for cid in ids:
        match = pattern.match(cid)
        if match:
            plate_name = match.group(1) or "Plate_1"
            row_letter = match.group(2).upper()
            col_number = int(match.group(3))
            well_id = f"{row_letter}{col_number:02d}"
        else:
            plate_name = "Plate_1"
            well_id = cid
            row_letter = np.nan
            col_number = np.nan

        plates.append(plate_name)
        wells.append(well_id)
        rows.append(row_letter)
        cols.append(col_number)

    if "plate" not in adata.obs:
        adata.obs["plate"] = plates
    adata.obs["well"] = wells
    adata.obs["well_row"] = rows
    adata.obs["well_col"] = cols

    console.print(
        f"[bold green]✔ Plate metadata parsed:[/bold green] {len(set(plates))} plate(s) detected."
    )
    return adata


def filter_cells(
    adata: ad.AnnData,
    min_genes: int = 1500,
    max_genes: Optional[int] = None,
    min_counts: int = 50000,
    max_counts: Optional[int] = None,
    max_pct_mito: float = 15.0,
    max_pct_ercc: Optional[float] = 15.0,
    max_doublet_score: Optional[float] = None,
) -> ad.AnnData:
    """Filter SMART-seq2 cells based on read depth, gene count, %mito, and %ERCC spike-ins."""
    n_initial = adata.n_obs
    mask = (adata.obs["n_genes_by_counts"] >= min_genes) & (adata.obs["total_counts"] >= min_counts)

    if max_genes is not None:
        mask &= adata.obs["n_genes_by_counts"] <= max_genes

    if max_counts is not None:
        mask &= adata.obs["total_counts"] <= max_counts

    if "pct_counts_mito" in adata.obs:
        mask &= adata.obs["pct_counts_mito"] <= max_pct_mito

    if max_pct_ercc is not None and "pct_counts_ercc" in adata.obs:
        mask &= adata.obs["pct_counts_ercc"] <= max_pct_ercc

    if max_doublet_score is not None and "doublet_score" in adata.obs:
        mask &= adata.obs["doublet_score"] <= max_doublet_score

    filtered_adata = adata[mask].copy()
    n_filtered = filtered_adata.n_obs
    console.print(
        f"[bold green]SMART-seq2 Cell Filtering:[/bold green] Retained {n_filtered}/{n_initial} cells "
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
    expected_doublet_rate: float = 0.02,
    random_state: int = 42,
    n_prin_comps: Optional[int] = None,
) -> ad.AnnData:
    """Run Scrublet doublet detection (optional for SMART-seq2 FACS-sorted plates)."""
    import scrublet as scr

    n_samples, n_vars = adata.shape
    if n_prin_comps is not None:
        n_pcs = min(n_prin_comps, n_samples - 1, n_vars - 1)
    else:
        n_pcs = min(30, max(2, (min(n_samples, n_vars) // 6)))

    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    console.print(
        f"[bold cyan]Running Scrublet doublet prediction (rate: {expected_doublet_rate}, pcs: {n_pcs})...[/bold cyan]"
    )

    scrub = scr.Scrublet(
        counts, expected_doublet_rate=expected_doublet_rate, random_state=random_state
    )
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=1,
        min_cells=1,
        min_gene_variability_pctl=85,
        n_prin_comps=n_pcs,
        verbose=False,
    )

    if predicted_doublets is None:
        predicted_doublets = np.zeros(adata.n_obs, dtype=bool)

    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets

    console.print(
        f"[bold green]✔ Scrublet finished:[/bold green] {predicted_doublets.sum()}/{adata.n_obs} predicted doublets."
    )
    return adata
