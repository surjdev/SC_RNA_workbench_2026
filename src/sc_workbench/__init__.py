"""
Single-Cell Downstream Analysis Workbench (sc-workbench).
Production-ready Python toolkit for downstream transcriptomics analysis.
"""

from .workbench import SingleCellWorkbench
from . import io
from . import qc
from . import preprocess
from . import reduction
from . import clustering
from . import markers
from . import annotation
from . import pathway
from . import trajectory
from . import plotting

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
