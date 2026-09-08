# %% [markdown]
# # Downstream analysis: AnnData + Pandas + NumPy + Scikit-learn
# รัน **Restart Kernel and Run All** ได้ด้วยข้อมูลสังเคราะห์ที่กำหนด seed ไม่มีการดาวน์โหลดข้อมูล
# เปลี่ยน `COUNTS_PATH` เพื่อใช้ matrix จาก Bash หรือ Nextflow (แถว = genes, คอลัมน์ = cells)
# ตัวอย่างนี้สาธิตการทำงานของซอฟต์แวร์ ไม่ใช่ข้อมูลสำหรับสรุปชีววิทยา

# %%
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

# Works from repository root, workbench/, or workbench/notebooks/.
root = next(
    p
    for p in [Path.cwd(), *Path.cwd().parents]
    if (p / "sc_workbench").is_dir() or (p / "workbench/sc_workbench").is_dir()
)
package_root = root if (root / "sc_workbench").is_dir() else root / "workbench"
sys.path.insert(0, str(package_root))
from sc_workbench import SingleCellWorkbench, clustering, io, markers, preprocess, qc, reduction

sc.settings.verbosity = 0
GTF_PATH = None  # Set to the exact upstream genes.gtf when using Ensembl IDs
COUNTS_PATH = None  # e.g. Path('/absolute/path/results/counts/gene_cell_count_matrix.tsv')
OUTPUT = package_root / "results" / "tutorial"
OUTPUT.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# ## 1. Input และ metadata
# Pandas/NumPy ใช้ cells × genes; upstream TSV ใช้ genes × cells และ loader จะ transpose ให้
# ชื่อเซลล์ต้องไม่ซ้ำ ใช้ index ในการจับคู่ metadata เสมอ
# สำหรับ Ensembl IDs ให้ใส่ gene symbols จาก GTF รุ่นเดียวกันใน `adata.var['gene_name']` ก่อนคำนวณ mitochondrial QC

# %%
if COUNTS_PATH is None:
    rng = np.random.default_rng(42)
    truth = np.repeat(["A", "B", "C"], 40)
    rates = np.full((120, 300), 1.5)
    for i in range(3):
        rates[i * 40 : (i + 1) * 40, 10 + i * 30 : 40 + i * 30] = 9
    df = pd.DataFrame(
        rng.poisson(rates),
        index=[f"cell_{i:03}" for i in range(120)],
        columns=["MT-DEMO"] + [f"gene_{i}" for i in range(299)],
    )
    adata = io.from_dataframe(df)
    adata.obs["synthetic_group"] = truth
else:
    adata = io.load_upstream_matrix(COUNTS_PATH)
if GTF_PATH is not None:
    io.annotate_from_gtf(adata, GTF_PATH)
    print("Matched annotation:", adata.var["annotation_matched"].mean())
# Example real metadata: adata.obs = metadata.loc[adata.obs_names].copy()
adata

# %% [markdown]
# ## 2. QC และ filtering
# ตรวจ distribution แล้วเลือก threshold ตาม protocol ตัวเลขนี้มีไว้สำหรับ demo เท่านั้น
# ไม่เรียก Scrublet อัตโนมัติสำหรับ plate-based Smart-seq2; พิจารณาวิธีตรวจ doublet ตามการเตรียมตัวอย่าง
# QC ใช้ `layers['counts']` แม้ `.X` ผ่าน normalization แล้ว

# %%
qc.calculate_qc_metrics(adata)
display(adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mito"]].describe())
adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mito"]].hist(bins=20, figsize=(10, 3))
adata = qc.filter_cells(adata, min_genes=20, min_counts=50, max_pct_mito=30)
qc.filter_genes(adata, min_cells=3)
assert adata.n_obs > 3 and adata.n_vars > 3, "Inspect QC thresholds: too few cells/genes remain"

# %% [markdown]
# ## 3. Normalization → HVG → PCA → neighbors → UMAP → Leiden
# ฟังก์ชันส่วนใหญ่แก้ AnnData เดิมและคืน object เดิม; `filter_cells` คืนสำเนาที่กรองแล้ว
# `counts` = original counts, `normalized` และ `.raw` = log1p normalized expression
# เก็บ genes ทั้งหมดไว้สำหรับ markers; PCA เลือกเฉพาะ HVG ไม่จำเป็นต้อง scale ทุก gene จน dense
# [Scanpy HVG documentation](https://scanpy.readthedocs.io/en/stable/api/scanpy.pp.highly_variable_genes.html): `seurat` ใช้ log expression ส่วน `seurat_v3` ใช้ counts และต้องติดตั้ง scikit-misc เพิ่ม

# %%
preprocess.normalize_and_log(adata)
preprocess.select_hvg(adata, n_top_genes=min(150, adata.n_vars), flavor="seurat")
reduction.run_pca(adata, n_comps=20)
reduction.compute_neighbors(adata, n_neighbors=10, n_pcs=20, use_rep="X_pca")
reduction.run_umap(adata, random_state=42)
clustering.cluster_leiden(adata, resolution=0.5)
sc.pl.umap(adata, color=["leiden", "total_counts"], show=True)

# %% [markdown]
# ## 4. ใช้ Pandas / NumPy / Scikit-learn โดยตรง
# ไม่จำเป็นต้องห่อ estimator ทุกตัวใน workbench ใช้ embedding ที่เก็บใน `.obsm` ได้ทันที
# ใช้ sparse matrix ผ่าน `.layers` สำหรับ matrix ใหญ่; `.to_numpy()` / `.to_df()` จะแปลงเป็น dense
# ARI ด้านล่างเทียบกับ label สังเคราะห์เท่านั้น ไม่ใช่ validation ของ cell type จริง
# หากฝึก classifier จริง ให้แยก train/test ตาม donor และ fit preprocessing เฉพาะ train เพื่อป้องกัน leakage

# %%
pca = io.get_embedding_df(adata, "X_pca")
adata.obs["kmeans"] = pd.Categorical(
    KMeans(n_clusters=3, random_state=42, n_init=10).fit_predict(pca).astype(str)
)
display(pd.crosstab(adata.obs["leiden"], adata.obs["kmeans"]))
if "synthetic_group" in adata.obs:
    print("Synthetic ARI:", adjusted_rand_score(adata.obs["synthetic_group"], adata.obs["kmeans"]))
# Inspect a small dense slice only
expression = io.to_dataframe(adata[:, :5], layer="normalized")
display(expression.head())
print("Counts per cell:", np.asarray(adata.layers["counts"].sum(axis=1)).ravel()[:5])

# %% [markdown]
# ## 5. Marker exploration และ checkpoint
# ใช้ log normalized expression สำหรับ marker tests ไม่ใช้ scaled matrix
# cluster markers เป็น exploratory; การเปรียบเทียบ treatment ต้องคำนึงถึง biological replicates
# อย่าตั้งชื่อ cell type จากหมายเลข cluster เพียงอย่างเดียว ต้องตรวจ marker และบริบท tissue

# %%
if adata.obs["leiden"].nunique() > 1:
    markers.find_markers(adata, groupby="leiden", n_genes=10)
    marker_table = markers.get_markers_df(adata)
    display(marker_table.groupby("cluster", observed=True).head(3))
    marker_table.to_csv(OUTPUT / "markers.csv", index=False)
io.save_h5ad(adata, OUTPUT / "processed.h5ad")
adata.obs.to_csv(OUTPUT / "cell_metadata.csv")
pca.to_csv(OUTPUT / "pca.csv")
loaded = sc.read_h5ad(OUTPUT / "processed.h5ad")
assert loaded.shape == adata.shape
assert "counts" in loaded.layers

# %% [markdown]
# ## 6. Optional fluent API
# เป็น wrapper บาง ๆ บน AnnData เท่านั้น ใช้ `wb.adata` ร่วมกับ Scanpy ได้
# constructor รับ object เดิม หากต้องการแยกการแก้ไขให้ `.copy()` ก่อน

# %%
wb = SingleCellWorkbench(adata.copy())
wb.add_obs(pd.Series("reviewed", index=wb.obs.index[::-1]), col_name="review_status")
display(wb.get_embedding("X_umap").head())
wb
