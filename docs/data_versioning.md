# Data Version Control & Large File Management Workflow

คู่มือการจัดการเวอร์ชันของข้อมูลและไฟล์ขนาดใหญ่ในโปรเจค Single-Cell Transcriptomics Analysis Workbench ตามข้อกำหนด **FR-4** และ **NFR-6**

---

## 1. วัตถุประสงค์ (Objective)

ข้อมูล scRNA-seq เช่น raw count matrices, FASTQ, BAM, และไฟล์ AnnData (`.h5ad`), Loom (`.loom`) มักมีขนาดตั้งแต่หลายสิบ Megabytes ไปจนถึงหลายสิบ Gigabytes หาก commit เข้า Git โดยตรง จะทำให้ขนาดของ `.git` repository บวม ช้า และไม่สามารถ push ขึ้น GitHub ปกติได้

โปรเจคนี้ใช้ **Git Large File Storage (Git LFS)** เป็นแกนหลักในการดักจับและจัดการไฟล์ขนาดใหญ่ outside git tree โดยอัตโนมัติ

---

## 2. ไฟล์ที่ถูกตรวจจับและควบคุมด้วย Git LFS

กำหนดไว้ในไฟล์ [.gitattributes](file:///home/surj/Workspace/SC_RNA_workbench_2026/.gitattributes):

- `*.h5ad` (AnnData container files)
- `*.loom` (Loom format)
- `*.h5`, `*.hdf5` (HDF5 hierarchical datasets)
- `*.mtx.gz`, `*.tsv.gz`, `*.csv.gz` (Compressed count matrices)
- `*.bam`, `*.fastq.gz` (Raw sequencing alignments)

---

## 3. ขั้นตอนการทำงานสำหรับนักวิเคราะห์ (Analyst Workflow)

### 3.1 การเริ่มต้นใช้งานครั้งแรก (Initial Setup)

```bash
# ตรวจสอบและติดตั้ง Git LFS บนเครื่อง
git lfs install
```

### 3.2 การดึงข้อมูลขนาดใหญ่ลงมาทำงาน (Pulling Large Files)

เมื่อ clone repository ใหม่ ไฟล์ขนาดใหญ่จะถูกแทนที่ด้วย pointer file ขนาดเล็ก (~130 bytes) หากต้องการดึงไฟล์จริงลงมาในเครื่อง:

```bash
# ดึงไฟล์ LFS ทั้งหมด
git lfs pull

# หรือระบุดึงเฉพาะไฟล์ที่ต้องการ
git lfs pull --include="data/raw/pbmc3k.h5ad"
```

### 3.3 การเพิ่มไฟล์ข้อมูลใหม่ (Adding New Datasets)

วางไฟล์ลงใน [data/raw/](file:///home/surj/Workspace/SC_RNA_workbench_2026/data/raw/) หรือ [data/processed/](file:///home/surj/Workspace/SC_RNA_workbench_2026/data/processed/) แล้วสั่ง commit ตามปกติ Git LFS จะดักจับไฟล์อัตโนมัติ:

```bash
# 1. วางไฟล์ข้อมูล
cp /path/to/my_data.h5ad data/raw/my_data.h5ad

# 2. ตรวจสอบว่าไฟล์ถูกแทร็กผ่าน LFS หรือไม่
git check-attr -a data/raw/my_data.h5ad
# ผลลัพธ์ควรแสดง: filter: lfs

# 3. Stage และ commit pointer เข้า git
git add data/raw/my_data.h5ad
git commit -m "data: add PBMC 3k raw dataset"

# 4. Push ข้อมูลไปยัง remote repository
git push origin main
```

---

## 4. โครงสร้างและข้อตกลงการจัดเก็บข้อมูล (Data Folder Conventions)

| ไดเรกทอรี | หน้าที่ | นโยบายการจัดเก็บ |
|---|---|---|
| `data/raw/` | ข้อมูลดิบที่ได้รับมา (Raw matrices, 10x folders) | Read-only ห้ามแก้ไขไฟล์เดิมเด็ดขาด |
| `data/processed/` | ผลลัพธ์หลังผ่าน QC, Normalization, Batch correction (`.h5ad`) | สร้างและบันทึกอัตโนมัติจาก Notebooks |
| `configs/` | พารามิเตอร์และ metadata ที่ใช้ประมวลผล | Commit เข้า Git เป็น text file ปกติ |

---

## 5. ทางเลือกเสริม: DVC (Data Version Control with Cloud Storage)

ในกรณีที่ทีมต้องการเชื่อมต่อ Remote Storage ภายนอก (เช่น AWS S3, Google Cloud Storage, MinIO):
1. สามารถติดตั้ง DVC: `pixi add dvc` หรือ `pip install dvc`
2. เริ่มต้น DVC: `dvc init`
3. กำหนด remote storage:
   ```bash
   dvc remote add -d myremote s3://my-scrna-bucket/data
   dvc add data/raw/dataset.h5ad
   dvc push
   ```
