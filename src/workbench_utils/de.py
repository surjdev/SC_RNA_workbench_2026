"""
Differential Expression Module supporting Hybrid R/Python execution (FR-6, FR-11, FR-12, FR-13).
Provides isolated R interop layer (DESeq2, limma-voom via rpy2) alongside pure-Python PyDESeq2.
Complies with scRNAseq_Workbench_Requirements_v2.md.
"""

from typing import List, Optional, Tuple, Union

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from rich.console import Console

console = Console()


class RDependencyError(RuntimeError):
    """Raised with clear installation guidance when R or required Bioconductor packages are missing."""

    pass


def check_r_dependencies(
    required_packages: Optional[List[str]] = None,
) -> Tuple[bool, Optional[str]]:
    """Check if rpy2 and required R packages are installed and functional without crashing.

    Parameters
    ----------
    required_packages : list of str, optional
        List of R package names to verify (default is ['DESeq2']).

    Returns
    -------
    is_available : bool
        True if all dependencies are available, False otherwise.
    error_message : str or None
        Actionable installation guidance if unavailable.
    """
    if required_packages is None:
        required_packages = ["DESeq2"]

    try:
        import rpy2.robjects as robjects
    except (ImportError, ModuleNotFoundError) as e:
        msg = (
            f"rpy2 is not installed or could not be loaded: {e}\n"
            "To install R interop dependencies, run:\n"
            "  - Pixi: pixi install -e r-stats\n"
            "  - Conda/Mamba: mamba env update -f environment.yml\n"
            "  - Pip: pip install -e .[r-stats]"
        )
        return False, msg
    except Exception as e:
        msg = f"R initialization error via rpy2: {e}"
        return False, msg

    missing_packages = []
    for pkg in required_packages:
        try:
            res = robjects.r(f"requireNamespace('{pkg}', quietly = TRUE)")[0]
            if not res:
                missing_packages.append(pkg)
        except Exception:
            missing_packages.append(pkg)

    if missing_packages:
        pkgs_str = ", ".join(f"'{p}'" for p in missing_packages)
        conda_pkgs = " ".join(f"bioconductor-{p.lower()}" for p in missing_packages)
        msg = (
            f"Required R/Bioconductor package(s) missing: {missing_packages}\n"
            "To install the missing packages, choose one of:\n"
            "  1. Pixi environment: pixi install -e r-stats\n"
            f"  2. Conda/Mamba: mamba install -c bioconda {conda_pkgs}\n"
            f"  3. Within R: if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(c({pkgs_str}))"
        )
        return False, msg

    return True, None


def prepare_pseudobulk_data(
    adata: ad.AnnData,
    design_factor: str,
    layer: Optional[str] = "counts",
    min_cells_per_gene: int = 3,
    sample_key: Optional[str] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Extract and format raw integer counts and sample metadata for pseudobulk DE analysis.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData containing single-cell counts.
    design_factor : str
        Column in adata.obs representing the primary biological comparison variable.
    layer : str, optional
        Layer to extract raw counts from (default is 'counts'). If not found, uses .X.
    min_cells_per_gene : int, default 3
        Pre-filter genes expressed in fewer than this number of cells/samples.
    sample_key : str, optional
        Column in adata.obs to aggregate cells into pseudo-bulk samples per biological replicate.
        If None, each cell is treated as an individual observation.

    Returns
    -------
    counts_df : pd.DataFrame
        Samples as rows, Genes as columns (integer values).
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
            factor_val = adata.obs.loc[idx, design_factor].iloc[0]
            pb_meta.append({sample_key: sid, design_factor: factor_val})

        counts_df = pd.DataFrame(pb_counts, index=unique_samples, columns=adata.var_names)
        clinical_df = pd.DataFrame(pb_meta, index=unique_samples)
    else:
        counts_df = pd.DataFrame(X_int, index=adata.obs_names, columns=adata.var_names)
        clinical_df = adata.obs[[design_factor]].copy()

    # 3. Filter genes with low expression
    gene_mask = (counts_df > 0).sum(axis=0) >= min_cells_per_gene
    counts_df = counts_df.loc[:, gene_mask]

    console.print(
        f"[bold green]✔ Pseudobulk data prepared:[/bold green] {counts_df.shape[0]} samples, "
        f"{counts_df.shape[1]} genes (design: '{design_factor}')."
    )
    return counts_df, clinical_df


# Backward-compatible alias
prepare_pydeseq2_data = prepare_pseudobulk_data


def deseq2_r(
    counts_df: pd.DataFrame,
    clinical_df: pd.DataFrame,
    design_factors: Union[str, List[str]] = "condition",
    contrast: Optional[Tuple[str, str, str]] = None,
    quiet: bool = True,
    r_seed: int = 42,
) -> pd.DataFrame:
    """Execute differential expression using R's reference DESeq2 implementation via rpy2 (FR-11).

    Parameters
    ----------
    counts_df : pd.DataFrame
        Counts matrix (Samples as rows and Genes as columns, or Genes x Samples).
    clinical_df : pd.DataFrame
        Metadata dataframe indexed matching sample IDs.
    design_factors : str or list of str
        Factor(s) for design formula (e.g. "condition" -> "~ condition").
    contrast : tuple of (str, str, str), optional
        Contrast specification: (factor_name, test_level, reference_level).
    quiet : bool, default True
        Suppress intermediate DESeq2 messages.
    r_seed : int, default 42
        Random seed for R reproducibility (NFR-1).

    Returns
    -------
    pd.DataFrame
        Tidy results containing gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj.
    """
    ok, err_msg = check_r_dependencies(["DESeq2"])
    if not ok:
        raise RDependencyError(err_msg)

    import rpy2.robjects as robjects
    from rpy2.robjects import conversion, pandas2ri
    from rpy2.robjects.conversion import localconverter

    # Ensure countData in R is (Genes x Samples)
    if list(counts_df.index) == list(clinical_df.index):
        counts_matrix = counts_df.T
    elif list(counts_df.columns) == list(clinical_df.index):
        counts_matrix = counts_df.copy()
    else:
        common_samples = [
            s for s in clinical_df.index if s in counts_df.index or s in counts_df.columns
        ]
        if not common_samples:
            raise ValueError("No matching sample IDs between counts_df and clinical_df index!")
        if set(common_samples).issubset(counts_df.index):
            counts_matrix = counts_df.loc[common_samples].T
            clinical_df = clinical_df.loc[common_samples]
        else:
            counts_matrix = counts_df[common_samples]
            clinical_df = clinical_df.loc[common_samples]

    counts_matrix = counts_matrix.round().astype(int)

    if isinstance(design_factors, str):
        factors = [design_factors]
    else:
        factors = list(design_factors)
    formula_str = f"~ {' + '.join(factors)}"

    console.print(
        f"[bold cyan]Executing R DESeq2 via rpy2 (formula: {formula_str}, seed: {r_seed})...[/bold cyan]"
    )

    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_counts = conversion.py2rpy(counts_matrix)
        r_coldata = conversion.py2rpy(clinical_df)

    robjects.globalenv["r_counts_matrix"] = r_counts
    robjects.globalenv["r_clinical_df"] = r_coldata
    robjects.r(f"set.seed({r_seed})")

    r_script = f"""
    suppressPackageStartupMessages(library(DESeq2))
    dds <- DESeqDataSetFromMatrix(
        countData = r_counts_matrix,
        colData = r_clinical_df,
        design = as.formula("{formula_str}")
    )
    dds <- DESeq(dds, quiet = {"TRUE" if quiet else "FALSE"})
    """

    if contrast is not None:
        c_factor, c_test, c_ref = contrast
        r_script += f"""
        res <- results(dds, contrast = c("{c_factor}", "{c_test}", "{c_ref}"))
        """
    else:
        r_script += """
        res <- results(dds)
        """

    r_script += """
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df
    """

    r_res = robjects.r(r_script)

    with localconverter(robjects.default_converter + pandas2ri.converter):
        results_df = conversion.rpy2py(r_res)

    results_df = results_df.rename(
        columns={
            "log2FoldChange": "log2FoldChange",
            "lfcSE": "lfcSE",
            "stat": "stat",
            "pvalue": "pvalue",
            "padj": "padj",
            "baseMean": "baseMean",
        }
    )

    std_cols = ["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    ordered_cols = [c for c in std_cols if c in results_df.columns] + [
        c for c in results_df.columns if c not in std_cols
    ]
    results_df = results_df[ordered_cols]

    if "padj" in results_df.columns:
        results_df = results_df.sort_values(by="padj", ascending=True).reset_index(drop=True)

    console.print(
        f"[bold green]✔ R DESeq2 completed:[/bold green] {len(results_df)} genes analyzed."
    )
    return results_df


def limma_voom_r(
    counts_df: pd.DataFrame,
    clinical_df: pd.DataFrame,
    design_factor: str = "condition",
    contrast: Optional[Tuple[str, str, str]] = None,
    quiet: bool = True,
) -> pd.DataFrame:
    """Execute differential expression using R's limma-voom reference implementation via rpy2.

    Parameters
    ----------
    counts_df : pd.DataFrame
        Counts matrix (Samples as rows and Genes as columns, or Genes x Samples).
    clinical_df : pd.DataFrame
        Metadata dataframe indexed matching sample IDs.
    design_factor : str, default 'condition'
        Factor to test.
    contrast : tuple of (str, str, str), optional
        Contrast specification: (factor_name, test_level, reference_level).
    quiet : bool, default True
        Suppress intermediate messages.

    Returns
    -------
    pd.DataFrame
        Tidy results containing gene, baseMean, log2FoldChange, stat, pvalue, padj.
    """
    ok, err_msg = check_r_dependencies(["limma", "edgeR"])
    if not ok:
        raise RDependencyError(err_msg)

    import rpy2.robjects as robjects
    from rpy2.robjects import conversion, pandas2ri
    from rpy2.robjects.conversion import localconverter

    if list(counts_df.index) == list(clinical_df.index):
        counts_matrix = counts_df.T
    else:
        counts_matrix = counts_df.copy()

    counts_matrix = counts_matrix.round().astype(int)

    console.print(
        f"[bold cyan]Executing R limma-voom via rpy2 (factor: {design_factor})...[/bold cyan]"
    )

    with localconverter(robjects.default_converter + pandas2ri.converter):
        r_counts = conversion.py2rpy(counts_matrix)
        r_coldata = conversion.py2rpy(clinical_df)

    robjects.globalenv["r_counts_limma"] = r_counts
    robjects.globalenv["r_coldata_limma"] = r_coldata

    r_script = f"""
    suppressPackageStartupMessages({{
        library(edgeR)
        library(limma)
    }})
    dge <- DGEList(counts = r_counts_limma)
    dge <- calcNormFactors(dge)
    group <- as.factor(r_coldata_limma${design_factor})
    design <- model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    v <- voom(dge, design, plot = FALSE)
    fit <- lmFit(v, design)
    """

    if contrast is not None:
        _, test_lvl, ref_lvl = contrast
        r_script += f"""
        contr <- makeContrasts({test_lvl} - {ref_lvl}, levels = design)
        fit2 <- contrasts.fit(fit, contr)
        fit2 <- eBayes(fit2)
        top_res <- topTable(fit2, number = Inf, sort.by = "P")
        """
    else:
        r_script += """
        fit2 <- eBayes(fit)
        top_res <- topTable(fit2, coef = ncol(design), number = Inf, sort.by = "P")
        """

    r_script += """
    top_df <- as.data.frame(top_res)
    top_df$gene <- rownames(top_df)
    top_df
    """

    r_res = robjects.r(r_script)

    with localconverter(robjects.default_converter + pandas2ri.converter):
        results_df = conversion.rpy2py(r_res)

    results_df = results_df.rename(
        columns={
            "logFC": "log2FoldChange",
            "AveExpr": "baseMean",
            "t": "stat",
            "P.Value": "pvalue",
            "adj.P.Val": "padj",
        }
    )

    std_cols = ["gene", "baseMean", "log2FoldChange", "stat", "pvalue", "padj"]
    ordered_cols = [c for c in std_cols if c in results_df.columns] + [
        c for c in results_df.columns if c not in std_cols
    ]
    results_df = results_df[ordered_cols].reset_index(drop=True)

    console.print(
        f"[bold green]✔ R limma-voom completed:[/bold green] {len(results_df)} genes analyzed."
    )
    return results_df


def deseq2_python(
    counts_df: pd.DataFrame,
    clinical_df: pd.DataFrame,
    design_factors: Union[str, List[str]] = "condition",
    contrast: Optional[Tuple[str, str, str]] = None,
    quiet: bool = True,
) -> pd.DataFrame:
    """Execute PyDESeq2 pipeline as an explicit pure-Python alternative to the R engine (FR-11).

    Parameters
    ----------
    counts_df : pd.DataFrame
        Counts matrix (Samples as rows and Genes as columns).
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

    if "padj" in results_df.columns:
        results_df = results_df.sort_values(by="padj", ascending=True)

    console.print(
        f"[bold green]✔ PyDESeq2 complete:[/bold green] {len(results_df)} genes analyzed."
    )
    return results_df


# Backward-compatible alias
run_pydeseq2 = deseq2_python


def run_de(
    counts_df: pd.DataFrame,
    clinical_df: pd.DataFrame,
    design_factors: Union[str, List[str]] = "condition",
    contrast: Optional[Tuple[str, str, str]] = None,
    engine: str = "deseq2_r",
    quiet: bool = True,
    r_seed: int = 42,
    **kwargs,
) -> pd.DataFrame:
    """Config-driven Differential Expression dispatcher supporting R DESeq2, PyDESeq2, and limma-voom.

    Parameters
    ----------
    counts_df : pd.DataFrame
        Counts matrix.
    clinical_df : pd.DataFrame
        Sample annotations.
    design_factors : str or list of str, default 'condition'
        Factor(s) for experimental design.
    contrast : tuple of (str, str, str), optional
        Contrast specification: (factor_name, test_level, reference_level).
    engine : str, default 'deseq2_r'
        DE engine choice:
        - 'deseq2_r': Native R DESeq2 via rpy2 (reference implementation)
        - 'pydeseq2': Pure-Python alternative implementation
        - 'limma_voom_r': Native R limma-voom via rpy2
    quiet : bool, default True
        Suppress intermediate engine logs.
    r_seed : int, default 42
        Seed for numerical reproducibility.

    Returns
    -------
    pd.DataFrame
        Standardized differential expression results.
    """
    engine_normalized = engine.lower().strip()

    if engine_normalized in ["deseq2_r", "r_deseq2", "r-deseq2", "deseq2"]:
        return deseq2_r(
            counts_df=counts_df,
            clinical_df=clinical_df,
            design_factors=design_factors,
            contrast=contrast,
            quiet=quiet,
            r_seed=r_seed,
        )
    elif engine_normalized in ["pydeseq2", "deseq2_python", "python_deseq2"]:
        return deseq2_python(
            counts_df=counts_df,
            clinical_df=clinical_df,
            design_factors=design_factors,
            contrast=contrast,
            quiet=quiet,
        )
    elif engine_normalized in ["limma", "limma_voom", "limma_voom_r"]:
        factor = design_factors if isinstance(design_factors, str) else design_factors[0]
        return limma_voom_r(
            counts_df=counts_df,
            clinical_df=clinical_df,
            design_factor=factor,
            contrast=contrast,
            quiet=quiet,
        )
    else:
        raise ValueError(
            f"Unsupported differential expression engine: '{engine}'. "
            "Supported options: 'deseq2_r' (default R engine), 'pydeseq2' (pure Python), 'limma_voom_r'."
        )
