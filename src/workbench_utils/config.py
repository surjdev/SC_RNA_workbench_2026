"""
Configuration Management and Validation Helper Module.
Loads, validates, and supplies sensible defaults for config-driven single-cell workflows.
Complies with scRNAseq_Workbench_Requirements.md (FR-5, NFR-1).
"""

import sys
from pathlib import Path
from typing import Any, Dict, Optional, Union

import yaml
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

console = Console()

DEFAULT_CONFIG: Dict[str, Any] = {
    "project": {
        "name": "SMARTseq2_Workbench_Analysis",
        "random_seed": 42,
    },
    "data": {
        "input_path": "data/raw/smartseq2_counts.tsv",
        "metadata_path": "data/raw/plate_metadata.csv",
        "output_h5ad": "data/processed/analyzed_workbench.h5ad",
        "report_dir": "reports",
    },
    "qc": {
        "min_genes": 1500,
        "max_genes": 12000,
        "min_counts": 50000,
        "max_counts": 10000000,
        "max_pct_mito": 15.0,
        "max_pct_ribo": 40.0,
        "max_pct_ercc": 15.0,
        "run_doublet_detection": False,
        "expected_doublet_rate": 0.02,
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
        "engine": "deseq2_r",  # Hybrid v2 default: R DESeq2 reference implementation via rpy2
        "tool": "deseq2_r",  # Alias for backward compatibility
        "r_seed": 42,
        "design_factor": "condition",
        "min_cells_per_gene": 3,
        "fdr_cutoff": 0.05,
        "log2fc_cutoff": 1.0,
    },
    "plotting": {
        "palette": "nature",
        "dpi": 150,
        "save_figures": True,
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
    config = {k: (v.copy() if isinstance(v, dict) else v) for k, v in DEFAULT_CONFIG.items()}

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

    validate_config(config)
    return config


def validate_config(config: Dict[str, Any]) -> bool:
    """Validate that required sections and critical types exist in configuration."""
    required_sections = ["project", "qc", "normalization", "reduction", "clustering"]
    for sec in required_sections:
        if sec not in config:
            raise KeyError(f"Missing required configuration section: '{sec}'")

    if "random_seed" not in config["project"]:
        raise KeyError("Missing 'random_seed' in config['project'] (required for NFR-1)!")

    if "differential_expression" in config:
        de_cfg = config["differential_expression"]
        engine = de_cfg.get("engine") or de_cfg.get("tool")
        valid_engines = {
            "deseq2_r",
            "r_deseq2",
            "deseq2",
            "pydeseq2",
            "deseq2_python",
            "limma_voom_r",
            "limma",
        }
        if engine and engine.lower().strip() not in valid_engines:
            raise ValueError(
                f"Invalid differential expression engine '{engine}'. "
                f"Must be one of: {sorted(list(valid_engines))}"
            )

    return True


def save_config(config: Dict[str, Any], output_path: Union[str, Path]) -> None:
    """Export configuration dictionary to a YAML file."""
    p = Path(output_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)
    console.print(f"[bold green]✔ Saved configuration to:[/bold green] {p}")


def print_config_summary(config: Dict[str, Any]) -> None:
    """Display pretty formatted table of configuration sections and values."""
    table = Table(title=f"Configuration: {config.get('project', {}).get('name', 'Pipeline')}")
    table.add_column("Section", style="cyan", no_wrap=True)
    table.add_column("Key", style="magenta")
    table.add_column("Value", style="green")

    for sec, values in config.items():
        if isinstance(values, dict):
            for k, v in values.items():
                table.add_row(sec, str(k), str(v))
        else:
            table.add_row(sec, "-", str(values))

    console.print(table)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        cfg_file = sys.argv[1]
        try:
            cfg = load_config(cfg_file)
            print_config_summary(cfg)
            console.print(Panel("[bold green]Validation Successful![/bold green]"))
        except Exception as e:
            console.print(Panel(f"[bold red]Validation Error:[/bold red] {e}"))
            sys.exit(1)
    else:
        print_config_summary(DEFAULT_CONFIG)
