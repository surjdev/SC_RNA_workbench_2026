"""
Configuration Management and Validation Helper Module.
Loads, validates, and supplies sensible defaults for config-driven single-cell workflows.
Complies with scRNAseq_Workbench_Requirements.md (FR-5, NFR-1).
"""

from pathlib import Path
from typing import Any, Dict, Optional, Union

import yaml
from rich.console import Console

console = Console()

DEFAULT_CONFIG: Dict[str, Any] = {
    "project": {
        "name": "scRNAseq_Workbench_Analysis",
        "random_seed": 42,
    },
    "qc": {
        "min_genes": 200,
        "max_genes": 8000,
        "min_counts": 500,
        "max_counts": 60000,
        "max_pct_mito": 20.0,
        "max_pct_ribo": 40.0,
        "expected_doublet_rate": 0.06,
        "max_doublet_score": 0.35,
    },
    "filter_genes": {
        "min_cells": 3,
    },
    "normalization": {
        "target_sum": 10000.0,
        "log1p": True,
        "n_top_genes": 2000,
        "flavor": "seurat",
    },
    "reduction": {
        "n_pcs": 30,
        "n_neighbors": 15,
        "metric": "cosine",
    },
    "clustering": {
        "method": "leiden",
        "resolution": 0.8,
    },
    "differential_expression": {
        "design_factor": "condition",
        "min_cells_per_gene": 3,
        "fdr_cutoff": 0.05,
        "log2fc_cutoff": 1.0,
    },
}


def load_config(config_path: Optional[Union[str, Path]] = None) -> Dict[str, Any]:
    """Load analysis YAML configuration file, merging user overrides with default parameters.

    Parameters
    ----------
    config_path : str or Path, optional
        Path to YAML configuration file. If None, returns default configuration.

    Returns
    -------
    dict
        Merged configuration dictionary.
    """
    config = DEFAULT_CONFIG.copy()

    if config_path is not None:
        p = Path(config_path)
        if not p.exists():
            raise FileNotFoundError(f"Configuration file not found: {p.resolve()}")

        with open(p, "r", encoding="utf-8") as f:
            user_config = yaml.safe_load(f) or {}

        # Deep merge section by section
        for section, values in user_config.items():
            if isinstance(values, dict) and section in config:
                config[section] = {**config[section], **values}
            else:
                config[section] = values

        console.print(f"[bold green]✔ Configuration loaded from:[/bold green] {p.name}")

    return config
