"""workbench_utils: Lightweight, non-invasive utility module for SMART-seq2 transcriptomics.

Provides shared I/O, QC helpers with ERCC spike-ins, plate layout plotting, config loaders, and PyDESeq2 utilities.
Complies with scRNAseq_Workbench_Requirements.md (FR-6, NFR-3, NFR-4, NFR-5).
"""

from workbench_utils.config import DEFAULT_CONFIG, load_config
from workbench_utils.de import (
    RDependencyError,
    check_r_dependencies,
    deseq2_python,
    deseq2_r,
    limma_voom_r,
    prepare_pseudobulk_data,
    prepare_pydeseq2_data,
    run_de,
    run_pydeseq2,
)
from workbench_utils.io import (
    export_for_seurat,
    load_h5ad,
    load_smartseq2_matrix,
    load_upstream_matrix,
    save_h5ad,
)
from workbench_utils.plotting import (
    CELL_PALETTE,
    NATURE_PALETTE,
    plot_plate_layout,
    plot_qc_violins,
    plot_volcano,
    set_publication_style,
)
from workbench_utils.qc import (
    calculate_qc_metrics,
    detect_doublets,
    filter_cells,
    filter_genes,
    parse_plate_metadata,
)

__version__ = "0.2.0"

__all__ = [
    "DEFAULT_CONFIG",
    "load_config",
    "prepare_pydeseq2_data",
    "prepare_pseudobulk_data",
    "run_pydeseq2",
    "deseq2_python",
    "deseq2_r",
    "limma_voom_r",
    "run_de",
    "check_r_dependencies",
    "RDependencyError",
    "load_upstream_matrix",
    "load_smartseq2_matrix",
    "load_h5ad",
    "save_h5ad",
    "export_for_seurat",
    "set_publication_style",
    "plot_qc_violins",
    "plot_plate_layout",
    "plot_volcano",
    "NATURE_PALETTE",
    "CELL_PALETTE",
    "calculate_qc_metrics",
    "parse_plate_metadata",
    "filter_cells",
    "filter_genes",
    "detect_doublets",
]
