# sc_workbench: downstream analysis in Python

เริ่มที่ [notebooks/00_end_to_end.ipynb](notebooks/00_end_to_end.ipynb): notebook ครบในไฟล์เดียว
พร้อม synthetic data, QC, normalization, HVG, PCA, UMAP, Leiden, Scikit-learn, markers และ H5AD export
notebooks 01–05 เป็นตัวอย่างเดิมสำหรับหัวข้อเพิ่มเติม ไม่ได้เป็นเส้นทาง smoke test หลัก

## Installation

```bash
cd workbench
pixi run lab
# หรือ virtual environment
python -m pip install -e '.[notebooks,dev]'
python -m jupyterlab
```

สำหรับ script นอก directory นี้ ให้ติดตั้ง `pip install -e /path/to/workbench` ใน Python environment ที่ใช้
Jupyter kernel ต้องใช้ environment เดียวกับที่ติดตั้ง dependencies

## Data model

| ตำแหน่ง | ความหมาย |
|---|---|
| `adata.X` | cells × genes; เปลี่ยนตาม processing step |
| `adata.layers['counts']` | input expression ที่เก็บตอน ingest; สำหรับงาน count pipeline ให้ใช้ raw counts เท่านั้น |
| `adata.layers['normalized']` | library normalized + log1p expression |
| `adata.raw` | snapshot log expression เมื่อ normalize(save_raw=True); ไม่ใช่ raw counts |
| `adata.obs` / `adata.var` | Pandas cell / gene metadata |
| `adata.obsm['X_pca']` / `['X_umap']` | embeddings สำหรับ plotting / Scikit-learn |

Pandas/NumPy constructors ใช้ cells × genes; upstream TSV loader ใช้ genes × cells โดย default
ชื่อเซลล์และ gene IDs ต้อง unique; metadata จับคู่ตาม index ไม่ใช่ลำดับแถว
ใช้ `io.annotate_from_gtf(adata, 'genes.gtf')` เติม gene_name/chromosome โดยจับคู่ gene_id แบบตรงตัว
ตรวจ `adata.var['annotation_matched']` ว่าจับคู่ได้ครบก่อนวิเคราะห์
QC ใช้ counts layer และ gene_name สำหรับตรวจ MT/RPS/RPL prefixes ถ้ามี
หาก input เป็น Ensembl IDs ต้องเติม symbols จาก annotation รุ่นเดียวกับ reference ก่อน mitochondrial filtering

```python
from sc_workbench import io, qc, preprocess, reduction
from sklearn.cluster import KMeans
import pandas as pd

adata = io.load_upstream_matrix('gene_cell_count_matrix.tsv')
qc.calculate_qc_metrics(adata)
# Inspect adata.obs distributions before choosing thresholds.
adata = qc.filter_cells(adata, min_genes=200, min_counts=500, max_pct_mito=20)
qc.filter_genes(adata, min_cells=3)
preprocess.normalize_and_log(adata)
preprocess.select_hvg(adata, n_top_genes=2000)
reduction.run_pca(adata, n_comps=30)
adata.obs['kmeans'] = pd.Categorical(
    KMeans(n_clusters=3, random_state=42, n_init=10)
    .fit_predict(adata.obsm['X_pca']).astype(str)
)
io.save_h5ad(adata, 'results/processed.h5ad')
```

ฟังก์ชันส่วนใหญ่ mutate AnnData และคืน object เดิม; `qc.filter_cells` คืน filtered copy
ต้องการทดลองแยก branch ให้ใช้ `adata.copy()`
`SingleCellWorkbench(adata)` เป็น optional fluent wrapper ที่อ้าง object เดิม ใช้ `wb.adata` ร่วมกับ Scanpy ได้ตลอด

`wb.add_obs(series, 'column')` align ตาม cell IDs; `.to_df()` / `.to_numpy()` คืน dense matrix
สำหรับข้อมูลใหญ่ ใช้ `adata.layers['counts']` เป็น sparse โดยตรง หรือ slice genes ก่อน export
การระบุ layer/representation ผิดจะ raise error ไม่มี fallback ไปข้อมูลอื่น
`normalize()` เริ่มจาก counts layer ทุกครั้ง จึงไม่ normalize log expression ซ้ำ
`seurat_v3` HVG ใช้ counts และต้องติดตั้ง `scikit-misc` เพิ่ม; default `seurat` ใช้ log expression

## ขอบเขตการตีความ

QC thresholds, resolution, batch correction และ annotation ต้องเลือกตาม dataset
Scrublet เป็น optional และจะ raise error ถ้าทำงานไม่ได้ ไม่มีการสร้างผลว่า “ไม่เป็น doublet” แทน failure
`train_classifier` คืน cross-validation report แต่ preprocessing ที่ทำทั้ง dataset ก่อน CV ยังมี leakage ได้
สำหรับ prediction จริงใช้ sklearn Pipeline และ split ตาม biological replicate/donor;
ความแม่นยำทำนาย Leiden labels ไม่ใช่หลักฐานยืนยัน cell types
pathway enrichment ต้องเลือก organism, identifiers และ background ให้ตรง; Enrichr ต้องใช้ network
Pseudotime ต้องกำหนด root จากชีววิทยา ไม่ใช่หมายเลข cluster โดยอัตโนมัติ

```bash
pixi run test
```

หลักการ HVG: [Scanpy documentation](https://scanpy.readthedocs.io/en/stable/api/scanpy.pp.highly_variable_genes.html)
