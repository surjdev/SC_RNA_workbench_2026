# %% [markdown]
# Start with [00_end_to_end.ipynb](00_end_to_end.ipynb) for an executable synthetic walkthrough. This topic notebook requires the indicated input checkpoints and dataset-specific parameters. Open from workbench/notebooks/.

# %% [markdown]
# # 🌿 04. Pathway Enrichment & Lineage Trajectory Inference
#
# This notebook demonstrates Over-Representation Analysis (Enrichr via gseapy), PAGA connectivity graph abstraction, and Diffusion Pseudotime (DPT) trajectory inference.

# %%
import sys
from pathlib import Path

sys.path.insert(0, str(Path("..").resolve()))

import scanpy as sc

import sc_workbench as scw
from sc_workbench import SingleCellWorkbench

wb = SingleCellWorkbench.from_h5ad("../data/annotated.h5ad")
wb

# %% [markdown]
# ### 1. Pathway Enrichment via Enrichr (GSEAPY)

# %%
# Run Enrichr on marker genes from cluster 0
markers_df = scw.markers.get_markers_df(wb.adata)
cluster0_genes = list(markers_df[markers_df["cluster"] == "0"]["gene"][:30])

# Requires real gene symbols, measured-gene background and network.
RUN_ENRICHR = False
if RUN_ENRICHR and len(cluster0_genes) > 0:
    enr_res = scw.pathway.run_enrichr(
        cluster0_genes,
        gene_sets="GO_Biological_Process_2023",
        organism="human",
        background=wb.adata.var_names.tolist(),
    )
    if not enr_res.empty:
        display(enr_res[["Term", "Overlap", "Adjusted P-value", "Genes"]].head(10))

# %% [markdown]
# ### 2. PAGA Lineage Graph Abstraction

# %%
# Set a root cell only after reviewing biological evidence.
ROOT_CELL = None
if ROOT_CELL is not None:
    scw.trajectory.run_paga(wb.adata, groups="leiden")
    sc.pl.paga(wb.adata, color="leiden", show=False)
    scw.trajectory.run_dpt(wb.adata, root_cell=ROOT_CELL)
    scw.plotting.plot_umap(wb.adata, color=["leiden", "dpt_pseudotime"])

# %% [markdown]
# ### 3. Diffusion Pseudotime (DPT) Inference

# %%
scw.trajectory.run_dpt(wb.adata)
scw.plotting.plot_umap(wb.adata, color=["cell_type", "dpt_pseudotime"])

# %% [markdown]
# ### 4. Export Seurat-Compatible Files

# %%
wb.export_seurat("../data/seurat_export")
print("✔ Stage 04 completed successfully!")
