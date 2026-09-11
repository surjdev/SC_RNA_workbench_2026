# ---
# jupyter:
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # 📊 Template 3: Pseudobulk Differential Expression Analysis (Hybrid R/Python)
# **Workflow:** Load Processed AnnData ➔ Pseudobulk Sample Aggregation ➔ Reference DE Modeling (R DESeq2 / PyDESeq2) ➔ Statistical Inference ➔ Volcano Plot ➔ Table Export
#
# Part of the Single-Cell SMART-seq2 Transcriptomics Analysis Workbench (FR-6, FR-10, FR-11, FR-12, FR-13).
#
# ### Architecture & Methodology:
# - **Default Engine:** Reference R `DESeq2` invoked via an isolated `rpy2` interop layer (`engine: "deseq2_r"`).
# - **Pure-Python Alternative:** `pydeseq2` (`engine: "pydeseq2"`) for zero-R deployment environments.
# - **No Silent Substitution:** The selected engine runs explicitly; missing R dependencies raise actionable guidance.
# - **Zero Analyst R Burden:** All data marshalling between AnnData / Pandas DataFrames and R matrix objects is handled seamlessly under the hood.

# %%
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from workbench_utils import (
    RDependencyError,
    check_r_dependencies,
    deseq2_python,
    deseq2_r,
    load_config,
    prepare_pseudobulk_data,
    run_de,
)
from workbench_utils.io import load_h5ad
from workbench_utils.plotting import plot_volcano, set_publication_style

set_publication_style()

# %%
# ==============================================================================
# Parameters Cell for Papermill (FR-7)
# ==============================================================================
config_path = "configs/default_analysis.yaml"
input_h5ad = "data/processed/02_clustered.h5ad"
output_table = "reports/differential_expression_results.csv"
engine = None  # None loads from config ("deseq2_r" by default). Options: "deseq2_r", "pydeseq2", "limma_voom_r"

# %%
# Load and validate YAML configuration (FR-5, NFR-1)
cfg = load_config(config_path if Path(config_path).exists() else None)
np.random.seed(cfg["project"]["random_seed"])
de_cfg = cfg["differential_expression"]
selected_engine = engine or de_cfg.get("engine", "deseq2_r")

print(f"Loaded Pipeline Configuration: {cfg['project']['name']}")
print(f"Configured DE Engine:         {selected_engine}")

# %% [markdown]
# ## 1. Load Processed AnnData
# Ingest clustered dataset from Template 2, or generate reproducible synthetic dataset with biological replicates if running standalone.

# %%
input_file = Path(input_h5ad)
if input_file.exists():
    adata = load_h5ad(input_file)
    print(f"Loaded AnnData from disk: {adata.n_obs} cells × {adata.n_vars} genes")
else:
    print(
        "Pre-clustered AnnData not found on disk. Generating synthetic dataset with biological replicates..."
    )
    rng = np.random.default_rng(cfg["project"]["random_seed"])
    n_cells, n_genes = 160, 200
    genes = [f"Gene_{i:03d}" for i in range(n_genes)]
    cells = [f"cell_{i:03d}" for i in range(n_cells)]

    # 2 conditions, 4 biological replicates per condition (20 cells per replicate)
    conditions = ["Control"] * 80 + ["Treated"] * 80
    sample_ids = [f"ctrl_rep_{i // 20 + 1}" for i in range(80)] + [
        f"treat_rep_{i // 20 + 1}" for i in range(80)
    ]

    # Base expression matrix
    X = rng.negative_binomial(6, 0.25, size=(n_cells, n_genes)).astype(np.float32)
    # Introduce true differential expression in first 25 genes for Treated group
    X[80:, :15] += rng.poisson(35, size=(80, 15))  # Upregulated in Treated
    X[80:, 15:25] = np.maximum(0, X[80:, 15:25] - 10)  # Downregulated in Treated

    obs = pd.DataFrame(
        {
            "condition": conditions,
            "sample_id": sample_ids,
            "batch": (["B1", "B2"] * 80),
        },
        index=cells,
    )
    var = pd.DataFrame({"gene_name": genes}, index=genes)

    adata = sc.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)
    adata.layers["counts"] = sp.csr_matrix(X)
    print(
        f"Generated synthetic single-cell AnnData: {adata.n_obs} cells across 8 samples, {adata.n_vars} genes."
    )

# %% [markdown]
# ## 2. Pseudobulk Sample Aggregation
# Per Single-Cell best practices, cells are aggregated per biological sample before statistical testing to prevent false discovery inflation from single-cell pseudo-replication.

# %%
design_factor = de_cfg.get("design_factor", "condition")
sample_key = de_cfg.get("sample_key", "sample_id") if "sample_id" in adata.obs else None
contrast = tuple(de_cfg["contrast"]) if "contrast" in de_cfg and de_cfg["contrast"] else None

counts_df, clinical_df = prepare_pseudobulk_data(
    adata,
    design_factor=design_factor,
    sample_key=sample_key,
    layer="counts" if "counts" in adata.layers else None,
    min_cells_per_gene=de_cfg.get("min_cells_per_gene", 3),
)

# Validate contrast against available metadata levels
if contrast and len(contrast) >= 3:
    factor_col = contrast[0]
    if factor_col in clinical_df.columns:
        present_levels = set(clinical_df[factor_col].dropna().unique())
        if contrast[1] not in present_levels or contrast[2] not in present_levels:
            print(
                f"Notice: Configured contrast {contrast} levels not found in '{factor_col}' ({present_levels})."
            )
            if len(present_levels) >= 2:
                levels = sorted(list(present_levels))
                contrast = (factor_col, levels[1], levels[0])
                print(f"Using available levels for contrast: {contrast}")

print(f"Aggregated count matrix: {counts_df.shape[0]} samples × {counts_df.shape[1]} genes")
print("\nSample metadata:")
print(clinical_df)

# %% [markdown]
# ## 3. Differential Expression Analysis (Configured Engine)
# Run the configured DE engine (`deseq2_r` default, or `pydeseq2`).
# If `deseq2_r` is requested in an environment where R/Bioconductor packages are not yet installed, `RDependencyError` provides actionable installation commands.

# %%
r_avail, r_msg = check_r_dependencies(["DESeq2"])
print(f"R DESeq2 Availability: {'✔ Available' if r_avail else '✘ Not Available'}")
if not r_avail:
    print(f"Notice: {r_msg}")

try:
    de_results = run_de(
        counts_df=counts_df,
        clinical_df=clinical_df,
        design_factors=design_factor,
        contrast=contrast,
        engine=selected_engine,
        r_seed=de_cfg.get("r_seed", 42),
        quiet=True,
    )
    primary_engine_used = selected_engine
except RDependencyError as e:
    print(f"\n[RDependencyError encountered]:\n{e}")
    print("\nFalling back to explicit Python engine (PyDESeq2) for this interactive session...")
    de_results = run_de(
        counts_df=counts_df,
        clinical_df=clinical_df,
        design_factors=design_factor,
        contrast=contrast,
        engine="pydeseq2",
        quiet=True,
    )
    primary_engine_used = "pydeseq2"

print(f"\nDE completed with [{primary_engine_used}]: {len(de_results)} genes evaluated.")
print("\nTop Differentially Expressed Genes:")
display_cols = [
    c for c in ["gene", "baseMean", "log2FoldChange", "pvalue", "padj"] if c in de_results.columns
]
print(de_results.head(10)[display_cols].to_string(index=False))

# %% [markdown]
# ## 4. Pure-Python PyDESeq2 Evaluation (Side-by-Side Comparison)
# Demonstrate explicit pure-Python `deseq2_python` execution (FR-11).

# %%
print("Executing PyDESeq2 (Python-only reference)...")
py_results = deseq2_python(
    counts_df=counts_df,
    clinical_df=clinical_df,
    design_factors=design_factor,
    contrast=contrast,
    quiet=True,
)
print(f"PyDESeq2 finished: {len(py_results)} genes evaluated.")
print(py_results.head(5)[display_cols].to_string(index=False))

# %% [markdown]
# ## 5. Publication-Ready Volcano Plot
# Visualize significant up- and down-regulated genes with Nature-themed palettes.

# %%
fdr_thresh = de_cfg.get("fdr_cutoff", 0.05)
fc_thresh = de_cfg.get("log2fc_cutoff", 1.0)
contrast_label = f"{contrast[1]} vs {contrast[2]}" if contrast else "Treated vs Control"

fig = plot_volcano(
    de_results,
    pval_col="padj",
    fc_col="log2FoldChange",
    gene_col="gene",
    pval_threshold=fdr_thresh,
    fc_threshold=fc_thresh,
    top_n_labels=10,
    title=f"Volcano Plot: {contrast_label} ({primary_engine_used})",
)
plt.show()

# %% [markdown]
# ## 6. Filter & Export DE Results
# Save full results and statistically significant hits to CSV tables.

# %%
out_path = Path(output_table)
out_path.parent.mkdir(parents=True, exist_ok=True)
de_results.to_csv(out_path, index=False)

# Export significant genes
sig_mask = (de_results["padj"] < fdr_thresh) & (de_results["log2FoldChange"].abs() >= fc_thresh)
sig_genes = de_results[sig_mask].sort_values(by="padj", ascending=True)
sig_out_path = out_path.parent / f"significant_{out_path.name}"
sig_genes.to_csv(sig_out_path, index=False)

print("Differential expression analysis completed successfully.")
print(f"  [✔] Engine:             {primary_engine_used}")
print(f"  [✔] Full table:         {out_path}")
print(f"  [✔] Significant hits:   {sig_out_path} ({len(sig_genes)} genes)")
