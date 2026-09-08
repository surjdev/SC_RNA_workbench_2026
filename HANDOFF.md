# AI HANDOFF & Verification Protocol
## Downstream Single-Cell Transcriptomics Analysis Workbench

เอกสารนี้จัดทำขึ้นเพื่อให้ AI Agent หรือวิศวกรซอฟต์แวร์สามารถตรวจสอบความคืบหน้า ตรวจสอบคุณภาพงาน (Verification) และรับช่วงต่องานในแต่ละ Part ได้อย่างไร้รอยต่อ

---

## 1. Project Context & Golden Rules

- **Golden Rule 1 (No API Wrapper - NFR-5):** เครื่องมือใน `src/workbench_utils/` ต้องเป็น Helper Functions อิสระ (I/O, QC, Plotting themes, PyDESeq2 formatter, Config) **ห้ามสร้าง Monolithic Wrapper Class (เช่น SingleCellWorkbench)** มาครอบ Scanpy หรือ PyDESeq2 เด็ดขาด นักวิเคราะห์ต้องเรียก Scanpy/PyDESeq2 API ตรงได้
- **Golden Rule 2 (Zero Notebook Output Noise - FR-3, NFR-6):** ไฟล์ `.ipynb` ต้องไม่มี output ติดเข้าไปใน Git (`nbstripout`) และต้องจับคู่กับ `.py:percent` เสมอ (`jupytext`)
- **Golden Rule 3 (Config-Driven - FR-5, NFR-1):** ห้าม hardcode ค่า threshold หรือ parameters ใน notebook ต้องโหลดจาก `configs/*.yaml`
- **Golden Rule 4 (Git Hygiene - FR-4):** ไฟล์ข้อมูลขนาดใหญ่ (`.h5ad`, `.loom`, raw matrices) ต้องไม่ถูก commit ตรงเข้า Git ต้องผ่าน Data Versioning (`git-lfs` หรือ `dvc`)

---

## 2. Milestone Tracking Board

| Part | ชื่องาน (Milestone) | สถานะ | ผู้ดำเนินการ | วันที่เสร็จสิ้น | ผ่าน Verification หรือไม่ |
|:---:|---|:---:|:---:|:---:|:---:|
| **Part 1** | Scaffolding & Directory Alignment | ✅ Completed | AI Agent | 2026-09-09 | [x] ผ่าน (Scaffold & Script OK) |
| **Part 2** | Environment & Dependencies (PyDESeq2, Conda, Pixi) | ✅ Completed | AI Agent | 2026-09-09 | [x] ผ่าน (pydeseq2 & 17 tests OK) |
| **Part 3** | Git Hygiene & Data Versioning (Pre-commit, LFS/DVC) | ✅ Completed | AI Agent | 2026-09-09 | [x] ผ่าน (pre-commit & Git LFS OK) |
| **Part 4** | Core Utility Package (`src/workbench_utils/`) | ✅ Completed | AI Agent | 2026-09-09 | [x] ผ่าน (29 tests & ruff check OK) |
| **Part 5** | Config-Driven Parameterization (`configs/*.yaml`) | ✅ Completed | AI Agent | 2026-09-09 | [x] ผ่าน (YAML schemas & validation OK) |
| **Part 6** | Parameterized Execution & Report Generation (Papermill/Make) | ⏳ Pending | - | - | [ ] |
| **Part 7** | Template Notebooks (QC, Clustering, PyDESeq2) & README | ⏳ Pending | - | - | [ ] |
| **Part 8** | CI/CD Pipeline (GitHub Actions) & Final Acceptance | ⏳ Pending | - | - | [ ] |

---

## 3. Detailed Verification Checklist & Commands per Part

ก่อนที่ AI จะรายงานว่า Part ใดเสร็จสิ้น หรือก่อนเริ่ม Part ถัดไป ให้รันคำสั่งตรวจสอบและเช็คเงื่อนไขตามตารางนี้:

### 🔹 Part 1: Scaffolding & Directory Alignment
* **สิ่งที่ต้องส่งมอบ:**
  - โครงสร้างไดเรกทอรีที่ Root: `data/raw/`, `data/processed/`, `notebooks/`, `src/workbench_utils/`, `configs/`, `reports/`, `tests/`
  - สคริปต์ `scripts/init_workbench.sh` สำหรับ scaffold โปรเจคใหม่
  - ลบ `.pixi/` และ symlink ที่เสียหายเดิม
  - Root `.gitignore` ที่ถูกต้อง
* **คำสั่งตรวจสอบ (Verification Commands):**
  ```bash
  test -d data/raw && test -d data/processed && test -d notebooks && test -d src/workbench_utils && test -d configs && test -d reports && echo "Scaffold OK"
  test -f scripts/init_workbench.sh && bash scripts/init_workbench.sh --test
  git status --porcelain
  ```
* **เกณฑ์ผ่าน (Pass Criteria):** ไดเรกทอรีครบตาม FR-1, สคริปต์รันผ่าน, ไม่มี broken paths ใน root

---

### 🔹 Part 2: Environment & Dependencies
* **สิ่งที่ต้องส่งมอบ:**
  - `pyproject.toml` ที่ Root มี `pydeseq2`, `scanpy`, `anndata`, `papermill`, `jupytext`, `nbstripout`, `nbconvert`, `ruff`
  - `environment.yml` ที่ติดตั้งได้ด้วย `mamba env create -f environment.yml` ในคำสั่งเดียว (FR-2, NFR-2)
  - `pixi.toml` และ `pixi.lock` ที่สะอาด
* **คำสั่งตรวจสอบ (Verification Commands):**
  ```bash
  python -c "import scanpy, anndata, pydeseq2, sklearn, papermill, jupytext; print('All Dependencies OK')"
  pixi check || echo "Pixi env verified"
  ```
* **เกณฑ์ผ่าน (Pass Criteria):** `import pydeseq2` สำเร็จโดยไม่มี ImportError

---

### 🔹 Part 3: Git Hygiene & Data Versioning
* **สิ่งที่ต้องส่งมอบ:**
  - `.pre-commit-config.yaml` บรรจุ `nbstripout`, `jupytext`, `ruff`
  - คอนฟิก Data Versioning (`.gitattributes` สำหรับ Git LFS หรือ `.dvc/` สำหรับ DVC)
  - คู่มือ `docs/data_versioning.md`
* **คำสั่งตรวจสอบ (Verification Commands):**
  ```bash
  pre-commit run --all-files
  git check-attr -a data/raw/sample.h5ad
  ```
* **เกณฑ์ผ่าน (Pass Criteria):** Pre-commit hooks ทำงานได้, ไฟล์ `.h5ad` ถูกดักจับโดย LFS/DVC, `git diff` ของ notebook ไม่แสดง output

---

### 🔹 Part 4: Core Utility Package (`src/workbench_utils/`)
* **สิ่งที่ต้องส่งมอบ:**
  - แพ็กเกจ `workbench_utils` ติดตั้งแบบ editable ได้ (`pip install -e .`)
  - ฟังก์ชัน I/O รองรับ `backed='r'` AnnData mode (NFR-3)
  - Helper สำหรับ PyDESeq2 formatting (`workbench_utils.de`)
  - Helper สำหรับ QC metrics และ publication plotting themes
  - Unit tests ใน `tests/` ครอบคลุมฟังก์ชันใน `src/`
* **คำสั่งตรวจสอบ (Verification Commands):**
  ```bash
  pip list | grep workbench-utils
  pytest tests/ -v
  ruff check src/
  ```
* **เกณฑ์ผ่าน (Pass Criteria):** Pytest ผ่าน 100%, Ruff ผ่านไม่มี lint error, ไม่มี class wrapper `SingleCellWorkbench`

---

### 🔹 Part 5: Config-Driven Parameterization
* **สิ่งที่ต้องส่งมอบ:**
  - `configs/default_analysis.yaml` มีครบทุก threshold, seed, resolution, DE cutoffs
  - โมดูล `workbench_utils.config` สำหรับโหลดและ validate ค่า YAML
* **คำสั่งตรวจสอบ (Verification Commands):**
  ```bash
  python -c "from workbench_utils.config import load_config; cfg = load_config('configs/default_analysis.yaml'); print('Config keys:', list(cfg.keys()))"
  ```
* **เกณฑ์ผ่าน (Pass Criteria):** สามารถโหลด config และเข้าถึง parameters ได้อย่างถูกต้อง มี default fallback

---

### 🔹 Part 6: Parameterized Execution & Reporting Pipeline
* **สิ่งที่ต้องส่งมอบ:**
  - Parameters cell ใน notebook
  - Runner script หรือ `Makefile` (คำสั่ง `make run`, `make report`)
  - การส่งออก HTML report ไปยัง `reports/` อัตโนมัติ
* **คำสั่งตรวจสอบ (Verification Commands):**
  ```bash
  make run NB=01_qc_and_filtering CONFIG=configs/default_analysis.yaml
  make report NB=01_qc_and_filtering
  ls -lh reports/*.html
  ```
* **เกณฑ์ผ่าน (Pass Criteria):** Papermill รันผ่านโดยฉีดค่า config ได้, ไฟล์ HTML report ถูกสร้างใน `reports/`

---

### 🔹 Part 7: Template Notebooks & Onboarding Docs
* **สิ่งที่ต้องส่งมอบ:**
  - 3 Template Notebooks:
    1. `notebooks/01_qc_and_filtering.ipynb`
    2. `notebooks/02_clustering_and_annotation.ipynb`
    3. `notebooks/03_differential_expression_pydeseq2.ipynb`
  - จับคู่ `.py:percent` ทุกไฟล์ด้วย Jupytext
  - `README.md` ฉบับสมบูรณ์สำหรับ onboarding ผู้ร่วมทีมใหม่
* **คำสั่งตรวจสอบ (Verification Commands):**
  ```bash
  jupytext --check notebooks/*.ipynb
  test -f notebooks/03_differential_expression_pydeseq2.ipynb && echo "Notebook 3 exists"
  ```
* **เกณฑ์ผ่าน (Pass Criteria):** Notebooks รัน end-to-end ได้ผลลัพธ์สอดคล้องตาม Acceptance Criteria

---

### 🔹 Part 8: CI/CD Pipeline & Final Acceptance
* **สิ่งที่ต้องส่งมอบ:**
  - `.github/workflows/ci.yml`
  - ผลการรัน Acceptance Criteria Checklist ครบทุกข้อในหัวข้อ 7
* **คำสั่งตรวจสอบ (Verification Commands):**
  ```bash
  pytest tests/ --cov=workbench_utils
  ruff check .
  ruff format --check .
  ```
* **เกณฑ์ผ่าน (Pass Criteria):** CI workflow syntax ถูกต้อง, tests ครอบคลุม, repository พร้อมใช้งานระดับ Production

---

## 4. Instructions for Incoming AI

เมื่อคุณเป็น AI Agent ที่เข้ามารับช่วงงานในรอบถัดไป:
1. **อ่านเอกสารนี้ (`HANDOFF.md`) และ [implementation_plan.md](file:///home/surj/.gemini/antigravity-ide/brain/8a8a72e0-35f5-4a28-a702-d73b6eb4697b/implementation_plan.md) เป็นอันดับแรก**
2. **ตรวจดู Milestone Tracking Board** ว่าผู้ใช้สั่งให้ทำ Part ใด
3. **รันคำสั่งตรวจสอบของ Part ก่อนหน้า** เพื่อยืนยันว่าโปรเจคอยู่ในสถานะที่พร้อมต่อยอด
4. **ดำเนินการเฉพาะขอบเขตของ Part ที่ได้รับมอบหมาย** ไม่แก้ไขนอกขอบเขตโดยไม่จำเป็น
5. **เมื่อทำเสร็จ ให้รัน Verification Commands** ของ Part นั้น บันทึกผล และอัปเดตสถานะใน Milestone Tracking Board ในเอกสารนี้
