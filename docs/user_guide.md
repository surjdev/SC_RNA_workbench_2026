# SMART-seq2 Transcriptomics Workbench: Comprehensive User Guide
## คู่มือการติดตั้ง การตั้งค่าสภาพแวดล้อม และการนำเข้าข้อมูล SMART-seq2 เข้าสู่การวิเคราะห์จริง

เอกสารคู่มือฉบับนี้จัดทำขึ้นเพื่อแนะนำขั้นตอนการใช้งาน **SMART-seq2 Plate-Based Single-Cell Transcriptomics Analysis Workbench** ตั้งแต่ขั้นตอนแรกสุด (Installation) จนถึงการนำไฟล์ข้อมูลดิบ Full-Length Read Counts และ Plate/Well Metadata เข้าสู่ Pipeline การวิเคราะห์เชิงลึก

---

## สารบัญ (Table of Contents)

1. [ปรัชญาการออกแบบและสถาปัตยกรรม (Architecture Philosophy)](#1-ปรัชญาการออกแบบและสถาปัตยกรรม)
2. [การติดตั้งและตั้งค่าสภาพแวดล้อม (Environment Setup)](#2-การติดตั้งและตั้งค่าสภาพแวดล้อม)
   - [ทางเลือกที่ 1: ติดตั้งผ่าน Pixi (แนะนำ - เร็วและล็อกเวอร์ชัน 100%)](#21-ติดตั้งผ่าน-pixi-recommended)
   - [ทางเลือกที่ 2: ติดตั้งผ่าน Conda / Mamba (มาตรฐานชีวสารสนเทศ)](#22-ติดตั้งผ่าน-conda--mamba)
   - [การตรวจสอบความพร้อมของระบบ (System Verification)](#23-การตรวจสอบความพร้อมของระบบ)
3. [คู่มือการนำเข้าข้อมูลดิบ SMART-seq2 (SMART-seq2 Data Ingestion Guide)](#3-คู่มือการนำเข้าข้อมูลดิบ-smart-seq2-smart-seq2-data-ingestion-guide)
   - [1. การนำเข้า Matrix ผลลัพธ์จาก Upstream Pipelines (STAR, RSEM, featureCounts, Kallisto)](#31-การนำเข้า-matrix-ผลลัพธ์จาก-upstream-pipelines)
   - [2. การนำเข้าข้อมูล Plate & Well Metadata (96/384-Well Layout)](#32-การนำเข้าข้อมูล-plate--well-metadata-96384-well-layout)
   - [3. การตรวจวัด Ambion ERCC Spike-in Controls และการแยกแยะหลุมที่ล้มเหลว (Failed Wells)](#33-การตรวจวัด-ambion-ercc-spike-in-controls-และการแยกแยะหลุมที่ล้มเหลว)
   - [4. การตรวจสอบรูปแบบเชิงพื้นที่ของ Plate (Plate Layout Heatmaps)](#34-การตรวจสอบรูปแบบเชิงพื้นที่ของ-plate-plate-layout-heatmaps)
   - [5. ไฟล์ AnnData (`.h5ad`) และ Loom (`.loom`)](#35-ไฟล์-anndata-h5ad-และ-loom-loom)
   - [6. การจัดการข้อมูลขนาดใหญ่ด้วย AnnData Backed Mode (`backed='r'`)](#36-การจัดการข้อมูลขนาดใหญ่ด้วย-anndata-backed-mode-backedr)
   - [7. การส่งออกข้อมูลข้ามภาษาไปยัง R / Seurat (`export_for_seurat`)](#37-การส่งออกข้อมูลข้ามภาษาไปยัง-r--seurat)
4. [การกำหนดค่าพารามิเตอร์ผ่านไฟล์ Config (YAML Configuration)](#4-การกำหนดค่าพารามิเตอร์ผ่านไฟล์-config)
5. [แนวทางการรันการวิเคราะห์ (Interactive vs Batch Report Execution)](#5-แนวทางการรันการวิเคราะห์)
   - [การทำงานแบบ Interactive ใน JupyterLab](#51-การทำงานแบบ-interactive-ใน-jupyterlab)
   - [การรันแบบอัตโนมัติด้วย Papermill และสร้าง HTML Report](#52-การรันแบบอัตโนมัติด้วย-papermill-และสร้าง-html-report)
6. [การควบคุมเวอร์ชันข้อมูลและ Git Hygiene](#6-การควบคุมเวอร์ชันข้อมูลและ-git-hygiene)
7. [คำถามที่พบบ่อยและการแก้ปัญหา (FAQ & Troubleshooting)](#7-คำถามที่พบบ่อยและการแก้ปัญหา)

---

## 1. ปรัชญาการออกแบบและสถาปัตยกรรม

Workbench นี้ถูกสร้างขึ้นเพื่อรองรับเทคโนโลยี **SMART-seq2 (Switching Mechanism at 5' End of RNA Template)** ซึ่งเป็นการจัดลำดับแบบ Plate-Based Full-Length Transcriptomics (96-well หรือ 384-well plates) โดยยึดหลักการสำคัญ:

- **No Custom Wrapper Class (NFR-5):** เครื่องมือใน `src/workbench_utils/` เป็น **Helper Functions อิสระ** ไม่มีคลาส Monolithic มาครอบ นักวิเคราะห์สามารถเรียกใช้ native APIs ของ `scanpy`, `pydeseq2`, `anndata`, `scikit-learn` ได้อย่างอิสระ 100%
- **Full-Length Deep Sequencing Biology:** รองรับความลึกการอ่าน 500,000 – 5,000,000 raw read counts ต่อเซลล์ และตรวจพบยีน 4,000 – 10,000 ยีนต่อเซลล์ แตกต่างจากเทคโนโลยี droplet ที่มีความลึกต่ำและ dropout สูง
- **Native PyDESeq2 Negative Binomial GLM:** เนื่องจาก SMART-seq2 เป็น Full-length read counts ข้อมูลจึงเข้ากันได้โดยตรงกับโมเดลการแจกแจงแบบ Negative Binomial ของ DESeq2 โดยไม่ต้องพึ่งพา pseudo-bulk approximation
- **Zero Notebook Output Noise (FR-3, NFR-6):** ทำงานร่วมกับ `jupytext` (จับคู่ `.ipynb` และ `.py:percent`) และ `nbstripout` ทำให้ Git diff มีเฉพาะบรรทัดโค้ดที่แก้ไข
- **Config-Driven Parameterization (FR-5, NFR-1):** พารามิเตอร์สำคัญทุกอย่าง (เช่น QC cutoffs, random seed, ERCC thresholds, clustering resolution, contrast) บรรจุอยู่ในไฟล์ `configs/*.yaml`

```text
SC_RNA_workbench_2026/
├── configs/                     # YAML configuration files (SMART-seq2 QC cutoffs, seeds, contrasts)
│   ├── default_analysis.yaml    # มาตรฐาน SMART-seq2 (CPM 1M, ERCC max 15%, counts > 50k)
│   ├── strict_qc_analysis.yaml  # High-stringency QC สำหรับ plate ที่มี dead cells สูง
│   └── smartseq2_plate_analysis.yaml # Multi-plate batch analysis
├── data/
│   ├── raw/                     # ไฟล์ข้อมูลดิบ (TSV, CSV, h5ad) Tracked by Git LFS
│   └── processed/               # ไฟล์ AnnData ที่ผ่าน QC, clustering, annotation
├── docs/                        # เอกสารคู่มือและ workflow documentation
├── notebooks/                   # Jupyter Notebooks จับคู่กับ Jupytext .py:percent
│   ├── 01_qc_and_filtering.py   # Plate ingestion, ERCC QC, plate layout heatmaps
│   ├── 02_clustering_and_annotation.py # CPM norm (exclude ERCC), PCA, UMAP, Leiden
│   └── 03_differential_expression_pydeseq2.py # PyDESeq2 GLM on raw full-length read counts
├── reports/                     # HTML reports สร้างอัตโนมัติจาก Papermill + figures
├── scripts/                     # Shell scripts (init_workbench.sh, run_notebook.sh, benchmark)
├── src/workbench_utils/         # Core helper functions (io, qc, de, plotting, config)
├── tests/                       # Unit tests & acceptance tests (pytest - 100% pass)
├── .gitattributes               # Git LFS routing rules
├── .pre-commit-config.yaml      # Automated code quality and notebook output stripping
├── environment.yml              # Conda/Mamba environment specification (≤ 2 commands)
├── pixi.toml / pixi.lock        # Pixi reproducible binary lockfile
├── pyproject.toml               # Python package metadata and tool configs
└── Makefile                     # Pipeline command shortcuts
```

---

## 2. การติดตั้งและตั้งค่าสภาพแวดล้อม

### 2.1 ติดตั้งผ่าน Pixi (Recommended)

Pixi เป็น package manager ยุคใหม่สำหรับวิทยาศาสตร์ข้อมูลที่รวดเร็วอย่างยิ่ง และล็อกเวอร์ชันระดับ binary lockfile (`pixi.lock`) ช่วยให้ผลลัพธ์เหมือนกันในทุกเครื่อง

1. **ติดตั้ง Pixi (หากยังไม่มีในเครื่อง):**
   ```bash
   curl -fsSL https://pixi.sh/install.sh | bash
   # เปิด terminal ใหม่ หรือ source ~/.bashrc
   ```

2. **Clone และเข้าสู่โฟลเดอร์โปรเจค:**
   ```bash
   cd /home/surj/Workspace/SC_RNA_workbench_2026
   ```

3. **ติดตั้ง Environment และ Dependencies (คำสั่งเดียว):**
   ```bash
   make setup
   # หรือรันด้วย pixi โดยตรง:
   # pixi install && pixi run install-editable && pixi run pre-commit install
   ```

4. **เปิดใช้งาน JupyterLab เพื่อเริ่มวิเคราะห์:**
   ```bash
   make lab
   # หรือ pixi run lab
   ```

---

### 2.2 ติดตั้งผ่าน Conda / Mamba

สำหรับทีมที่ใช้ระบบ Conda หรือ Mamba ในการจัดการ Environment:

1. **สร้าง Environment จาก `environment.yml` (คำสั่งเดียวตาม NFR-2):**
   ```bash
   mamba env create -f environment.yml
   # หรือ: conda env create -f environment.yml
   ```

2. **เปิดใช้งาน Environment:**
   ```bash
   conda activate sc_workbench
   ```

3. **ติดตั้ง workbench package ในโหมด editable:**
   ```bash
   pip install -e .
   pre-commit install
   ```

4. **ลงทะเบียน Jupyter Kernel:**
   ```bash
   python -m ipykernel install --user --name sc_workbench --display-name "Python 3.12 (sc_workbench)"
   ```

---

### 2.3 การตรวจสอบความพร้อมของระบบ

รันคำสั่งต่อไปนี้เพื่อยืนยันว่าการติดตั้งสมบูรณ์ 100%:

```bash
# ตรวจสอบการ import โมดูลสำคัญ
python -c "import scanpy, anndata, pydeseq2, workbench_utils; print('✔ Environment Ready!')"

# รัน unit tests ทั้งหมด
make test
# หรือ: pixi run pytest tests/ -v
```

---

## 3. คู่มือการนำเข้าข้อมูลดิบ SMART-seq2 (SMART-seq2 Data Ingestion Guide)

ฟังก์ชัน `workbench_utils.io.load_smartseq2_matrix` และ `load_upstream_matrix` ถูกออกแบบมาเพื่อนำเข้าผลลัพธ์จาก Upstream Processing Pipelines ของ SMART-seq2 ได้อย่างสมบูรณ์แบบ

### 3.1 การนำเข้า Matrix ผลลัพธ์จาก Upstream Pipelines

ในโปรโตคอล SMART-seq2 ขั้นตอน Upstream (Aligner + Quantifier) เช่น **STAR + RSEM**, **HISAT2 + featureCounts**, หรือ **Salmon / Kallisto** จะส่งออกเมทริกซ์การนับยีน (Gene Expression Count Matrix) ออกมาในรูปแบบ TSV หรือ CSV โดยปกติจะมีโครงสร้าง:
- **แถว (Rows):** รหัสยีนหรือสัญลักษณ์ยีน เช่น `ENSG00000...` หรือ `TP53`, รวมถึง Spike-in controls เช่น `ERCC-00002`
- **คอลัมน์ (Columns):** รหัสหลุมหรือเซลล์ เช่น `Plate1_A01`, `Plate1_A02`, ..., `Plate1_H12`

**วิธีนำเข้าด้วย `load_smartseq2_matrix`:**
```python
from workbench_utils.io import load_smartseq2_matrix

adata = load_smartseq2_matrix(
    count_path="data/raw/smartseq2_counts.tsv",
    metadata_path="data/raw/plate_metadata.csv",
    transpose=True,  # แปลงแถวยีนเป็นคอลัมน์ และคอลัมน์เซลล์เป็นแถว
    ercc_prefix=("ERCC-", "ercc-"),
)

print(adata)
# Output: AnnData object with n_obs × n_vars = 384 wells × 24,000 features
# adata.var['is_ercc'] จะเป็น True สำหรับยีนที่เป็น ERCC spike-in controls อัตโนมัติ
```

---

### 3.2 การนำเข้าข้อมูล Plate & Well Metadata (96/384-Well Layout)

ไฟล์ Plate Metadata (`plate_metadata.csv`) ช่วยระบุข้อมูลตำแหน่งเชิงกายภาพของหลุมในแผ่นทดสอบ ข้อมูลผู้บริจาค (Donor) และสภาวะการทดลอง (Condition):

| well_id | plate | well | well_row | well_col | condition | donor |
|---|---|---|---|---|---|---|
| Plate1_A01 | Plate1 | A01 | A | 1 | Treated | Donor_1 |
| Plate1_A02 | Plate1 | A02 | A | 2 | Treated | Donor_1 |
| Plate1_H12 | Plate1 | H12 | H | 12 | Control | Donor_2 |

หากไฟล์ Metadata มีคอลัมน์ `well` (เช่น `A01` ถึง `H12`) ระบบจะแยก `well_row` (ตัวอักษร A-H) และ `well_col` (ตัวเลข 1-12 หรือ 1-24) ให้โดยอัตโนมัติ เพื่อนำไปสร้าง Plate Layout Heatmap

---

### 3.3 การตรวจวัด Ambion ERCC Spike-in Controls และการแยกแยะหลุมที่ล้มเหลว

ใน SMART-seq2 จะมีการเติมสารควบคุมสังเคราะห์ **Ambion ERCC RNA Spike-In Control Mix** ปริมาณคงที่ลงในทุกหลุมก่อนทำ Reverse Transcription:
- **หลุมปกติ (Healthy single cells):** ปริมาณ endogenous cellular RNA มีสูงมาก ทำให้สัดส่วน ERCC spike-ins คิดเป็นเพียง **2% - 5%** ของ read ทั้งหมด
- **หลุมที่การแยกเซลล์ล้มเหลว (Empty Wells / FACS sorting failures):** ไม่มีเซลล์อยู่ในหลุม หรือเซลล์แตกสลายก่อน Reverse Transcription ทำให้ endogenous RNA ต่ำมาก ส่งผลให้สัดส่วน ERCC พุ่งสูงเกิน **20% - 80%** ของ read ทั้งหมด

```python
from workbench_utils.qc import calculate_qc_metrics, filter_cells

# 1. คำนวณ QC metrics รวมถึง pct_counts_ercc
calculate_qc_metrics(adata, ercc_prefix=("ERCC-", "ercc-"))

# 2. กรองหลุมที่ล้มเหลวทิ้ง (ค่า default: max_pct_ercc = 15.0%)
filtered_adata = filter_cells(adata, max_pct_ercc=15.0, min_counts=50000, min_genes=1500)
print(f"Retained wells: {filtered_adata.n_obs} / {adata.n_obs}")
```

---

### 3.4 การตรวจสอบรูปแบบเชิงพื้นที่ของ Plate (Plate Layout Heatmaps)

เพื่อตรวจสอบความผิดพลาดของหุ่นยนต์หยอดสาร (Liquid handler dispenser), การระเหยของของเหลวบริเวณขอบแผ่น (Edge effects), หรือปัญหาการเรียงลำดับหัว pipette:

```python
from workbench_utils.plotting import plot_plate_layout

# ตรวจสอบความลึกของการอ่านในแต่ละหลุมบน Plate 1 (Rows A-H, Cols 1-12)
fig_depth = plot_plate_layout(adata, plate_id="Plate1", color_key="total_counts")

# ตรวจสอบสัดส่วน ERCC spike-ins เพื่อดูตำแหน่งหลุมเปล่า
fig_ercc = plot_plate_layout(adata, plate_id="Plate1", color_key="pct_counts_ercc")
```

---

### 3.5 ไฟล์ AnnData (`.h5ad`) และ Loom (`.loom`)

หากได้รับไฟล์ที่เป็น AnnData หรือ Loom อยู่แล้ว:
```python
from workbench_utils.io import load_upstream_matrix

adata = load_upstream_matrix("data/raw/smartseq2_processed.h5ad", file_format="h5ad")
```

---

### 3.6 การจัดการข้อมูลขนาดใหญ่ด้วย AnnData Backed Mode (`backed='r'`)

เมื่อมีข้อมูลเซลล์ขนาดใหญ่มาก (เช่น > 200,000 เซลล์) ที่ไม่สามารถโหลดเข้า RAM ได้ทั้งหมดในคราวเดียว:
การเปิดไฟล์ด้วย `backed='r'` จะอ่านข้อมูลเฉพาะ metadata ลง RAM ส่วน count matrix ขนาดใหญ่จะยังคงอยู่บน Disk และจะถูกอ่านขึ้นมาเฉพาะตอนที่มีการ slice เท่านั้น (NFR-3):

```python
from workbench_utils.io import load_upstream_matrix

# เปิดไฟล์ h5ad ในโหมด read-only backed mode
adata_backed = load_upstream_matrix(
    path="data/raw/massive_atlas_500k_cells.h5ad",
    backed="r"
)

print("Is backed?", adata_backed.isbacked)  # True
print("Cell count:", adata_backed.n_obs)    # 500,000 cells ใช้ RAM เพียงไม่กี่ MB

# Slice เฉพาะเซลล์กลุ่มที่สนใจ (เช่น Subset เฉพาะ T Cells) แล้วโหลดเข้า Memory:
t_cell_mask = adata_backed.obs["cell_type"] == "T_Cell"
adata_t_cells = adata_backed[t_cell_mask].to_memory()
```

> **คำแนะนำ (Best Practice):** ใช้ `backed='r'` ในขั้นตอนเบื้องต้นเพื่อกรอง Sample หรือสำรวจ Metadata จากนั้นจึงใช้ `.to_memory()` โหลดเฉพาะ Subset ที่ต้องการวิเคราะห์เข้า RAM

---

### 3.7 การส่งออกข้อมูลข้ามภาษาไปยัง R / Seurat

หากต้องการส่งต่อข้อมูลที่ผ่านการ QC หรือ Clustering ใน Python ไปให้เพื่อนร่วมทีมวิเคราะห์ต่อใน R ด้วย Seurat:

**ใน Python:**
```python
from workbench_utils.io import export_for_seurat

# ส่งออกเป็น 10x MTX directory format + metadata.csv
export_for_seurat(
    adata=adata,
    output_dir="data/processed/seurat_export",
    layer="counts"  # หรือ None สำหรับ adata.X
)
```

**ใน R (Seurat):**
```R
library(Seurat)

# 1. อ่าน Matrix และ Metadata จาก CSV ที่ export จาก Python
counts <- read.csv("data/processed/seurat_export/counts.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("data/processed/seurat_export/metadata.csv", row.names = 1)

# 2. สร้าง Seurat Object
seurat_obj <- CreateSeuratObject(counts = as.matrix(counts), meta.data = meta, project = "SingleCellWorkbench")

# 3. หากมี embeddings.csv (UMAP)
if (file.exists("data/processed/seurat_export/embeddings.csv")) {
  umap_coords <- as.matrix(read.csv("data/processed/seurat_export/embeddings.csv", row.names = 1))
  seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = umap_coords, key = "UMAP_", assay = DefaultAssay(seurat_obj))
}

print(seurat_obj)
```

---

## 4. การกำหนดค่าพารามิเตอร์ผ่านไฟล์ Config

เพื่อหลีกเลี่ยงการ hardcode ค่า cutoffs และรักษา Reproducibility (NFR-1) ค่าพารามิเตอร์ทั้งหมดจะถูกจัดเก็บใน `configs/*.yaml`:

### ตัวอย่าง: `configs/default_analysis.yaml`
```yaml
project:
  name: "SMARTseq2_Default_Analysis"
  random_seed: 42

data:
  input_path: null
  metadata_path: null
  format: "smartseq2"
  transpose_counts: true

qc:
  min_counts: 50000        # SMART-seq2 read depth 500k-5M reads
  max_counts: 10000000
  min_genes: 1500          # ยีนที่ตรวจพบ > 1,500 - 4,000 ยีน
  max_genes: 12000
  mito_prefix: ["MT-", "mt-"]
  max_pct_mito: 15.0
  ribo_prefix: ["RPS", "RPL"]
  max_pct_ribo: 40.0
  ercc_prefix: ["ERCC-", "ercc-"]
  max_pct_ercc: 15.0       # กรองหลุมเปล่า/หลุมที่ FACS sort พลาด (ERCC > 15%)
  run_doublet_detection: false  # FACS sort 1 cell/well อัตรา doublet ต่ำมาก

normalization:
  target_sum: 1000000.0    # Counts Per Million (CPM)
  log1p: true
  n_top_genes: 2500
  flavor: "seurat"

reduction:
  n_pcs: 30
  n_neighbors: 15
  metric: "cosine"
  umap_min_dist: 0.5

clustering:
  resolution: 0.6
  algorithm: "leiden"

differential_expression:
  design_factor: "condition"
  contrast: ["condition", "Treated", "Control"]
  fdr_cutoff: 0.05
  log2fc_cutoff: 1.0
  min_cells_per_gene: 5
```

### การโหลดและ Validate Config ใน Python หรือ CLI:
```python
from workbench_utils.config import load_config

# โหลด config พร้อม fallback ค่า default อัตโนมัติ
cfg = load_config("configs/default_analysis.yaml")
seed = cfg["project"]["random_seed"]
min_genes = cfg["qc"]["min_genes"]
```

ตรวจสอบความถูกต้องของไฟล์ Config ผ่าน Command Line:
```bash
make validate-config CONFIG=configs/default_analysis.yaml
# หรือ: pixi run python -m workbench_utils.config configs/default_analysis.yaml
```

---

## 5. แนวทางการรันการวิเคราะห์

### 5.1 การทำงานแบบ Interactive ใน JupyterLab

เหมาะสำหรับการทดลองและสำรวจข้อมูลแบบ Real-time:
1. เปิด JupyterLab ด้วย `make lab`
2. เปิดไฟล์ในโฟลเดอร์ `notebooks/`:
   - `01_qc_and_filtering.ipynb`: SMART-seq2 Ingestion, ERCC Spike-ins & Mito QC Violins, Plate Layout Heatmaps, Failed Well Filtering
   - `02_clustering_and_annotation.ipynb`: Endogenous CPM Normalization (excluding ERCC), HVG, PCA, Plate Batch Inspection, UMAP, Leiden, Biomarkers
   - `03_differential_expression_pydeseq2.ipynb`: Raw full-length count extraction, PyDESeq2 Negative Binomial GLM Modeling, Volcano Plot, Table Export
3. เมื่อแก้ไขโค้ดใน JupyterLab ระบบ `jupytext` จะซิงค์การเปลี่ยนแปลงไปยังไฟล์ `.py` ที่เป็นคู่กันทันที

---

### 5.2 การรันแบบอัตโนมัติด้วย Papermill และสร้าง HTML Report

เหมาะสำหรับการรันซ้ำแบบ Batch Processing, การเปลี่ยนพารามิเตอร์ หรือการนำเข้าสู่ CI/CD:

```bash
# รัน Notebook 01 พร้อมฉีดค่า config และส่งออก HTML Report ไปยัง reports/
make report NB=01_qc_and_filtering CONFIG=configs/default_analysis.yaml

# รัน Notebook 02
make report NB=02_clustering_and_annotation CONFIG=configs/default_analysis.yaml

# รัน Notebook 03 สำหรับ Differential Expression
make report NB=03_differential_expression_pydeseq2 CONFIG=configs/default_analysis.yaml
```

รายงาน HTML จะถูกบันทึกใน `reports/<notebook_name>.html` ซึ่งสามารถเปิดดูผ่านเว็บเบราว์เซอร์เพื่อแชร์ให้ทีมนักวิจัยหรือบันทึกในสมุดบันทึกการทดลองได้ทันที

---

## 6. การควบคุมเวอร์ชันข้อมูลและ Git Hygiene

โปรเจคนี้ตั้งค่าระบบควบคุมความสะอาดของ Git เพื่อป้องกันไฟล์ขนาดใหญ่และ output ของ notebook รั่วไหลเข้าสู่ Git:

1. **การทำงานของ `nbstripout` และ `jupytext`:**
   - เมื่อสั่ง `git commit` ระบบ pre-commit จะลบ output ทั้งหมดออกจาก `.ipynb` อัตโนมัติ ทำให้ Git Diff แสดงเฉพาะโค้ดที่เปลี่ยนแปลงเท่านั้น
2. **การทำงานของ Git LFS:**
   - ไฟล์นามสกุล `.h5ad`, `.loom`, `.mtx.gz`, `.csv.gz` ใน `data/` จะถูกดักจับและส่งไปยัง Git LFS โดยอัตโนมัติตามกฎใน `.gitattributes`
   - Pre-commit hook `check-added-large-files` จะบล็อกไฟล์ขนาดเกิน 5 MB ไม่ให้ commit เข้า Git history ตรงๆ

---

## 7. คำถามที่พบบ่อยและการแก้ปัญหา (FAQ & Troubleshooting)

**Q1: ฉันจะเพิ่ม Notebook การวิเคราะห์ใหม่ได้อย่างไร?**
- สร้างไฟล์ใหม่ใน `notebooks/` เช่น `notebooks/04_trajectory_inference.ipynb`
- สั่งซิงค์ Jupytext เพื่อสร้างไฟล์ `.py` คู่กัน:
  ```bash
  pixi run jupytext --set-formats ipynb,py:percent notebooks/04_trajectory_inference.ipynb
  ```

**Q2: รัน `import pydeseq2` แล้วขึ้นแจ้งเตือน UserWarning เกี่ยวกับ `design_factors`?**
- ใน `workbench_utils.de.run_pydeseq2` มีการแปลง `design_factors` ให้เป็น formula string เช่น `~ condition` รองรับเวอร์ชัน 0.5.4+ โดยอัตโนมัติ จึงมั่นใจได้ว่าจะไม่เกิด deprecation error

**Q3: ฉันต้องการรันการทดสอบ Unit Tests ทั้งหมดในระบบ?**
- รันคำสั่งสั้นๆ:
  ```bash
  make test
  # หรือตรวจสอบ Code Coverage:
  pixi run pytest tests/ --cov=workbench_utils --cov-report=term-missing
  ```
