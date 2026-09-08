"""
Single-Cell Downstream Analysis Workbench (sc-workbench).
Production-ready Python toolkit for downstream transcriptomics analysis.
"""

from . import (
    annotation,
    clustering,
    io,
    markers,
    pathway,
    plotting,
    preprocess,
    qc,
    reduction,
    trajectory,
)
from .workbench import SingleCellWorkbench

__version__ = "0.1.0"
__all__ = [
    "SingleCellWorkbench",
    "io",
    "qc",
    "preprocess",
    "reduction",
    "clustering",
    "markers",
    "annotation",
    "pathway",
    "trajectory",
    "plotting",
]
