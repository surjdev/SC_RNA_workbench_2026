"""
workbench_utils: Lightweight, non-invasive utility module for single-cell transcriptomics.
Provides shared I/O, QC helpers, publication plotting themes, config loaders, and PyDESeq2 utilities.
Complies with scRNAseq_Workbench_Requirements.md (FR-6, NFR-3, NFR-4, NFR-5).
"""

from workbench_utils.config import DEFAULT_CONFIG, load_config
from workbench_utils.de import prepare_pydeseq2_data, run_pydeseq2
from workbench_utils.io import (
    export_for_seurat,
    load_10x_directory,
    load_h5ad,
    load_upstream_matrix,
    save_h5ad,
)
from workbench_utils.plotting import (
    CELL_PALETTE,
    NATURE_PALETTE,
    plot_qc_violins,
    plot_volcano,
    set_publication_style,
)
from workbench_utils.qc import (
    calculate_qc_metrics,
    detect_doublets,
    filter_cells,
    filter_genes,
)

__version__ = "0.1.0"

__all__ = [
    "DEFAULT_CONFIG",
    "load_config",
    "prepare_pydeseq2_data",
    "run_pydeseq2",
    "load_upstream_matrix",
    "load_10x_directory",
    "load_h5ad",
    "save_h5ad",
    "export_for_seurat",
    "set_publication_style",
    "plot_qc_violins",
    "plot_volcano",
    "NATURE_PALETTE",
    "CELL_PALETTE",
    "calculate_qc_metrics",
    "filter_cells",
    "filter_genes",
    "detect_doublets",
]
