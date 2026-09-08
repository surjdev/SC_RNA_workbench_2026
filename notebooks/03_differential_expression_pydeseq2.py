# %% [markdown]
# # 📊 Template 3: Differential Expression Analysis (PyDESeq2)
# **Workflow:** Load Clustered AnnData ➔ PyDESeq2 Model Fitting ➔ Statistical Inference ➔ Volcano Plot ➔ Table Export
#
# Part of the Single-Cell Transcriptomics Analysis Workbench (FR-10, NFR-5).
# Utilizes **PyDESeq2** Negative Binomial generalized linear models for rigorous DE testing.

# %%
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from workbench_utils.config import load_config
from workbench_utils.de import prepare_pydeseq2_data, run_pydeseq2
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

# %%
cfg = load_config(config_path if Path(config_path).exists() else None)
np.random.seed(cfg["project"]["random_seed"])
print(f"Loaded Pipeline Configuration: {cfg['project']['name']}")

# %% [markdown]
# ## 1. Load Processed AnnData
# Ingest clustered dataset from Template 2.

# %%
input_file = Path(input_h5ad)
if not input_file.exists():
    raise FileNotFoundError(
        f"Input file not found at: {input_file.resolve()}. Please run 02_clustering_and_annotation first!"
    )

adata = load_h5ad(input_file)
print(f"Loaded AnnData: {adata.n_obs} cells × {adata.n_vars} genes")

# %% [markdown]
# ## 2. Data Preparation for PyDESeq2
# Extract raw integer count matrix and sample metadata.
# Supports both single-cell testing and pseudo-bulk aggregation across biological replicates.

# %%
de_cfg = cfg["differential_expression"]
design_factor = de_cfg.get("design_factor", "condition")
sample_key = de_cfg.get("sample_key", None)

counts_df, clinical_df = prepare_pydeseq2_data(
    adata,
    design_factor=design_factor,
    layer="counts" if "counts" in adata.layers else None,
    min_cells_per_gene=de_cfg.get("min_cells_per_gene", 3),
    sample_key=sample_key,
)

print(f"Prepared count matrix: {counts_df.shape[0]} observations × {counts_df.shape[1]} genes")
print(f"Clinical metadata:\n{clinical_df.head()}")

# %% [markdown]
# ## 3. Fit PyDESeq2 Model & Statistical Testing
# Fit dispersion trend, run Negative Binomial GLM, and evaluate contrast.

# %%
contrast = tuple(de_cfg["contrast"]) if "contrast" in de_cfg and de_cfg["contrast"] else None

de_results = run_pydeseq2(
    counts_df=counts_df,
    clinical_df=clinical_df,
    design_factors=design_factor,
    contrast=contrast,
    quiet=True,
)

print(f"PyDESeq2 Testing complete: {len(de_results)} genes evaluated.")
print("\nTop Significant DE Genes:")
print(de_results.head(10)[["gene", "baseMean", "log2FoldChange", "pvalue", "padj"]])

# %% [markdown]
# ## 4. Publication-Ready Volcano Plot
# Visualize significant up- and down-regulated genes.

# %%
fdr_thresh = de_cfg.get("fdr_cutoff", 0.05)
fc_thresh = de_cfg.get("log2fc_cutoff", 1.0)

fig = plot_volcano(
    de_results,
    pval_col="padj",
    fc_col="log2FoldChange",
    gene_col="gene",
    pval_threshold=fdr_thresh,
    fc_threshold=fc_thresh,
    top_n_labels=10,
    title=f"PyDESeq2 Volcano Plot: {contrast[1] if contrast else 'Treated'} vs {contrast[2] if contrast else 'Control'}",
)
plt.show()

# %% [markdown]
# ## 5. Filter & Export DE Results
# Save full and significant DE gene tables to `reports/`.

# %%
out_path = Path(output_table)
out_path.parent.mkdir(parents=True, exist_ok=True)
de_results.to_csv(out_path, index=False)

# Export significant genes table
sig_mask = (de_results["padj"] < fdr_thresh) & (de_results["log2FoldChange"].abs() >= fc_thresh)
sig_genes = de_results[sig_mask]
sig_out_path = out_path.parent / f"significant_{out_path.name}"
sig_genes.to_csv(sig_out_path, index=False)

print("Differential expression analysis completed.")
print(f"  [✔] Full results saved to:        {out_path}")
print(f"  [✔] Significant genes ({len(sig_genes)}) saved to: {sig_out_path}")
