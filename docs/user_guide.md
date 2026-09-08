# Single-Cell Transcriptomics Workbench: Comprehensive User Guide
## คู่มือการติดตั้ง การตั้งค่าสภาพแวดล้อม และการนำเข้าข้อมูลสู่การวิเคราะห์จริง

เอกสารคู่มือฉบับนี้จัดทำขึ้นเพื่อแนะนำขั้นตอนการใช้งาน **Downstream Single-Cell Transcriptomics Analysis Workbench** ตั้งแต่ขั้นตอนแรกสุด (Installation) จนถึงการนำไฟล์ข้อมูลดิบทางชีววิทยาประเภทต่างๆ เข้าสู่ Pipeline การวิเคราะห์

---

## สารบัญ (Table of Contents)

1. [ปรัชญาการออกแบบและสถาปัตยกรรม (Architecture Philosophy)](#1-ปรัชญาการออกแบบและสถาปัตยกรรม)
2. [การติดตั้งและตั้งค่าสภาพแวดล้อม (Environment Setup)](#2-การติดตั้งและตั้งค่าสภาพแวดล้อม)
   - [ทางเลือกที่ 1: ติดตั้งผ่าน Pixi (แนะนำ - เร็วและล็อกเวอร์ชัน 100%)](#21-ติดตั้งผ่าน-pixi-recommended)
   - [ทางเลือกที่ 2: ติดตั้งผ่าน Conda / Mamba (มาตรฐานชีวสารสนเทศ)](#22-ติดตั้งผ่าน-conda--mamba)
   - [การตรวจสอบความพร้อมของระบบ (System Verification)](#23-การตรวจสอบความพร้อมของระบบ)
3. [คู่มือการนำเข้าข้อมูลดิบ (Data Ingestion Guide)](#3-คู่มือการนำเข้าข้อมูลดิบ-data-ingestion-guide)
   - [1. โฟลเดอร์ 10x Genomics Cell Ranger Matrix (`.mtx.gz`)](#31-โฟลเดอร์-10x-genomics-cell-ranger-matrix)
   - [2. ไฟล์ 10x Genomics HDF5 (`.h5`)](#32-ไฟล์-10x-genomics-hdf5-h5)
   - [3. ไฟล์ AnnData (`.h5ad`)](#33-ไฟล์-anndata-h5ad)
   - [4. ไฟล์ Loom Format (`.loom`)](#34-ไฟล์-loom-format-loom)
   - [5. ไฟล์ตารางข้อความ CSV / TSV Counts Matrix](#35-ไฟล์ตารางข้อความ-csv--tsv-counts-matrix)
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

Workbench นี้ถูกสร้างขึ้นภายใต้หลักการสำคัญตามเอกสารข้อกำหนด [scRNAseq_Workbench_Requirements.md](file:///home/surj/Workspace/SC_RNA_workbench_2026/scRNAseq_Workbench_Requirements.md):

- **No Custom Wrapper Class (NFR-5):** เครื่องมือใน `src/workbench_utils/` เป็น **Helper Functions อิสระ** ไม่มีคลาส Monolithic (เช่น `SingleCellWorkbench`) มาครอบ นักวิเคราะห์สามารถเรียกใช้ native APIs ของ `scanpy`, `pydeseq2`, `anndata`, `scikit-learn` ได้อย่างอิสระ 100%
- **Zero Notebook Output Noise (FR-3, NFR-6):** ทำงานร่วมกับ `jupytext` (จับคู่ `.ipynb` และ `.py:percent`) และ `nbstripout` ทำให้ Git diff มีเฉพาะบรรทัดโค้ดที่แก้ไข ไม่มี cell output ปะปน
- **Config-Driven Parameterization (FR-5, NFR-1):** พารามิเตอร์สำคัญทุกอย่าง (เช่น QC cutoffs, random seed, clustering resolution, contrast) ถูกแยกออกจากโค้ดและบรรจุอยู่ในไฟล์ `configs/*.yaml`

```text
SC_RNA_workbench_2026/
├── configs/                     # YAML configuration files (QC cutoffs, seeds, DE contrasts)
├── data/
│   ├── raw/                     # ไฟล์ข้อมูลดิบ (10x, h5ad, loom) Tracked by Git LFS
│   └── processed/               # ไฟล์ AnnData ที่ผ่าน QC, clustering, annotation
├── docs/                        # เอกสารคู่มือและ workflow documentation
├── notebooks/                   # Jupyter Notebooks จับคู่กับ Jupytext .py:percent
├── reports/                     # HTML reports สร้างอัตโนมัติจาก Papermill + figures
├── scripts/                     # Shell scripts (init_workbench.sh, run_notebook.sh)
├── src/workbench_utils/         # Core helper functions (io, qc, de, plotting, config)
├── tests/                       # Unit tests & acceptance tests (pytest)
├── .gitattributes               # Git LFS routing rules
├── .pre-commit-config.yaml      # Automated code quality and notebook output stripping
├── environment.yml              # Conda/Mamba environment specification (≤ 2 commands)
├── pixi.toml / pixi.lock        # Pixi reproducible dependency lockfile
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
   # หากไม่มี mamba สามารถใช้ conda env create -f environment.yml ได้เช่นกัน
   ```

2. **เปิดใช้งาน Environment:**
   ```bash
   conda activate sc_workbench
   ```

3. **ติดตั้งแพ็กเกจ `workbench_utils` แบบ Editable (`-e .`):**
   ```bash
   pip install -e .
   ```

4. **ติดตั้ง Git Pre-commit Hooks:**
   ```bash
   pre-commit install
   ```

5. **ลงทะเบียน Jupyter Kernel (เพื่อให้เลือก Kernel ใน Notebook ได้):**
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

## 3. คู่มือการนำเข้าข้อมูลดิบ (Data Ingestion Guide)

ฟังก์ชัน `workbench_utils.io.load_upstream_matrix` ถูกออกแบบมาเพื่อตรวจจับและโหลดไฟล์ข้อมูลดิบได้หลากหลายรูปแบบอัตโนมัติ โดยส่งคืนอ็อบเจกต์ `anndata.AnnData` มาตรฐานที่พร้อมส่งต่อให้ `scanpy` ใช้งานต่อได้ทันที

### 3.1 โฟลเดอร์ 10x Genomics Cell Ranger Matrix

ผลลัพธ์จาก Cell Ranger count มักจะอยู่ในโฟลเดอร์ `filtered_feature_bc_matrix/` ซึ่งประกอบด้วย 3 ไฟล์:
- `matrix.mtx.gz` (หรือ `matrix.mtx`)
- `barcodes.tsv.gz` (หรือ `barcodes.tsv`)
- `features.tsv.gz` (หรือ `genes.tsv.gz`)

**วิธีนำเข้า:**
```python
from workbench_utils.io import load_upstream_matrix

# ชี้ path ไปที่โฟลเดอร์ที่บรรจุ 3 ไฟล์ดังกล่าว
adata = load_upstream_matrix(
    path="data/raw/pbmc_sample/filtered_feature_bc_matrix",
    file_format="10x_mtx"
)

print(adata)
# Output: AnnData object with n_obs × n_vars = 2700 × 32738
```

หรือใช้ Scanpy Native API โดยตรง:
```python
import scanpy as sc

adata = sc.read_10x_mtx(
    "data/raw/pbmc_sample/filtered_feature_bc_matrix",
    var_names="gene_symbols",
    cache=True
)
```

---

### 3.2 ไฟล์ 10x Genomics HDF5 (`.h5`)

Cell Ranger ส่งออกไฟล์ `.h5` เช่น `filtered_feature_bc_matrix.h5` ซึ่งรวมทั้ง matrix และ metadata ไว้ในไฟล์เดียว:

**วิธีนำเข้า:**
```python
from workbench_utils.io import load_upstream_matrix

adata = load_upstream_matrix(
    path="data/raw/sample_filtered_feature_bc_matrix.h5",
    file_format="10x_h5"
)
```

หรือใช้ Scanpy Native API:
```python
import scanpy as sc

adata = sc.read_10x_h5("data/raw/sample_filtered_feature_bc_matrix.h5")
adata.var_names_make_unique()
```

---

### 3.3 ไฟล์ AnnData (`.h5ad`)

ไฟล์มาตรฐานของ Scanpy และ Single-Cell Python stack:

**วิธีนำเข้า:**
```python
from workbench_utils.io import load_upstream_matrix

adata = load_upstream_matrix(
    path="data/raw/reference_dataset.h5ad",
    file_format="h5ad"
)
```

หรือใช้ Scanpy Native API:
```python
import scanpy as sc

adata = sc.read_h5ad("data/raw/reference_dataset.h5ad")
```

---

### 3.4 ไฟล์ Loom Format (`.loom`)

ไฟล์ Loom มักพบในงานวิเคราะห์ RNA Velocity (velocyto) หรือข้อมูลจาก Human Cell Atlas:

**วิธีนำเข้า:**
```python
from workbench_utils.io import load_upstream_matrix

adata = load_upstream_matrix(
    path="data/raw/sample_velocity.loom",
    file_format="loom"
)
```

---

### 3.5 ไฟล์ตารางข้อความ CSV / TSV Counts Matrix

ในกรณีที่ได้รับข้อมูลนับดิบ (Raw Counts) จากงานวิจัยอื่นในรูปแบบ `.csv` หรือ `.tsv`:
- แถว (Rows): เซลล์ (Cell Barcodes) หรือ ยีน (Genes)
- คอลัมน์ (Columns): ยีน หรือ เซลล์

**วิธีนำเข้า:**
```python
from workbench_utils.io import load_upstream_matrix

adata = load_upstream_matrix(
    path="data/raw/counts_matrix.csv",
    file_format="csv"
)
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
  name: "scRNAseq_Default_Analysis"
  random_seed: 42

qc:
  min_counts: 500
  max_counts: 35000
  min_genes: 200
  max_genes: 6000
  max_pct_counts_mt: 20.0
  max_pct_counts_ribo: 50.0
  doublet_rate: 0.06

normalization:
  target_sum: 10000
  n_top_genes: 2000
  flavor: "seurat"

reduction:
  n_pcs: 30
  n_neighbors: 15

clustering:
  resolution: 0.8
  algorithm: "leiden"

differential_expression:
  design_factor: "condition"
  contrast: ["condition", "Treated", "Control"]
  padj_cutoff: 0.05
  log2fc_cutoff: 1.0
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
   - `01_qc_and_filtering.ipynb`: Data Ingestion, QC Violins, Scrublet Doublet Detection, Filtering
   - `02_clustering_and_annotation.ipynb`: Normalization, HVG, PCA, Harmony Integration, Leiden, UMAP, Cell Typing
   - `03_differential_expression_pydeseq2.ipynb`: Raw count extraction, PyDESeq2 GLM DE Analysis, Volcano Plot
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
