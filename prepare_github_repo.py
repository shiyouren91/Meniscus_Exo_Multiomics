"""
GitHub 仓库准备脚本 — 一键复制核心文件到 Meniscus_Exo_Multiomics/

使用方法:
  python prepare_github_repo.py

功能:
  1. 复制 R 脚本 (14个) → R/
  2. 复制 Python 核心脚本 → python/
  3. 复制 CSV 数据表 → data/tables/
  4. 复制投稿 PDF 图 → figures/
  5. 复制最终版投稿文档 → manuscript/
  6. 生成 DATA_BOM.md (物料清单)
  
注意: .rds 大文件不复制（通过 .gitignore 排除，或用 Git LFS）
"""

import os, shutil, hashlib, sys
from pathlib import Path
from datetime import datetime

# Force UTF-8 output on Windows
sys.stdout.reconfigure(encoding='utf-8')

# ============================================================
# 源路径配置
# ============================================================
SRC_R_ROOT = Path(r"C:\Users\89367\Desktop\CNS半月板")
SRC_PY_ROOT = Path(r"C:\Users\89367\WorkBuddy\20260331203952")
SRC_DATA_DIR = SRC_R_ROOT / "结果处理完" / "1ALL_RESULTS_FINAL"
SRC_SUPP_TABLES = Path(r"C:\Users\89367\Desktop\Nature_Submission_Meniscus_20260404\Supplementary_Tables")
SRC_FIGURES = Path(r"C:\Users\89367\Desktop\Nature_Submission_Meniscus_20260404\Figures")
SRC_MANUSCRIPT = Path(r"C:\Users\89367\Desktop\Nature_Submission_Meniscus_20260404\Manuscript")

# 目标仓库根目录
REPO_ROOT = Path(__file__).parent

# ============================================================
# 文件映射清单
# ============================================================

R_SCRIPTS = {
    # Module 1: scRNA integration
    "R/01_scRNA_integration": [
        ("meniscus_exosome_project/meniscus_exosome_project/scripts/01_scRNA_integration/step2_qc_integration.R", None),
        ("meniscus_exosome_project/meniscus_exosome_project/scripts/01_scRNA_integration/step3_annotation.R", None),
        ("meniscus_exosome_project/meniscus_exosome_project/scripts/01_scRNA_integration/step4_trajectory_GRN.R", None),
        ("meniscus_exosome_project/meniscus_exosome_project/scripts/01_scRNA_integration/step5_trajectory_degeneration.R", None),
    ],
    # Module 2: Exosome cargo
    "R/02_exosome_cargo": [
        ("meniscus_exosome_project/meniscus_exosome_project/scripts/02_exosome_cargo/step1_cargo_database.R", None),
        ("meniscus_exosome_project/meniscus_exosome_project/scripts/02_exosome_cargo/step2_miRNA_targets.R", None),
    ],
    # Module 3: Molecular docking
    "R/03_molecular_docking": [
        ("meniscus_exosome_project/meniscus_exosome_project/scripts/03_molecular_docking/step1_cargo_receptor_cellchat.R", None),
        ("meniscus_exosome_project/meniscus_exosome_project/scripts/03_molecular_docking/step2_cellchat_differential.R", None),
    ],
    # Module 4: Validation
    "R/04_validation": [
        ("meniscus_exosome_project/meniscus_exosome_project/scripts/04_validation/step1_cross_validation.R", None),
    ],
    # Auxiliary R scripts
    "R/auxiliary": [
        ("scenic_tf_analysis.R", None),           # Root level in CNS半月板
        ("trajectory_analysis.R", None),
        ("apply_annotation.R", None),
        ("module3_v2.R", None),
        ("validation_figs_v2.R", None),
    ],
}

PY_SCRIPTS = {
    "python/figures": [
        ("nc_figure_v5_antiartifact.py", "NC-standard figure rendering (Figures 1-4)"),
        ("nc_figure_v7_clean_rewrite.py", "Figure 4 v7 final rewrite"),
        ("publication_figures_script.py", "Batch publication figure generation"),
    ],
    "python/analysis": [
        ("trajectory_cargo_spatiotemporal.py", "Spatiotemporal cargo-trajectory mapping"),
        ("generate_supplementary.py", "Auto-generate Supplementary Information docx"),
    ],
    "python/submission_tools": [
        ("fix_maintext_v3_sync_suppmat.py", "Main Text v3 SI sync fixer"),
        ("nature_compliance_check.py", "NC format compliance checker"),
    ],
}

# CSV tables to copy from the canonical results directory
CSV_FILES = [
    "cluster_markers.csv",
    "cell_type_treatability_scores.csv",
    "protein_cargo_celltype_scores.csv",
    "miRNA_cargo_celltype_scores.csv",
    "TF_differential_activity.csv",
    "TF_exosome_cargo_links.csv",
    "TF_target_network.csv",
    "cargo_GO_enrichment.csv",
    "cargo_KEGG_enrichment.csv",
    "cargo_therapeutic_pathway_mapping.csv",
    "differentially_expressed_receptors.csv",
    "exosome_communication_repair_map.csv",
    "exosome_repair_map.csv",
    "hub_miRNAs.csv",
    "miRNA_high_confidence_targets.csv",
    "miRNA_target_network_edges.csv",
    "miRNA_targets_GO.csv",
    "MSC_source_cargo_comparison.csv",
    "optimal_cargo_recommendations.csv",
    "treatability_scores.csv",
    "druggability_assessment.csv",
    "validation_summary.csv",
    "analysis_summary.csv",
]

# Figures to copy
FIGURE_FILES = ["Figure1.pdf", "Figure1.png", "Figure2.pdf", "Figure2.png",
                "Figure3.pdf", "Figure3.png", "Figure4.pdf", "Figure4.png"]

# Manuscript files (final versions only)
MANUSCRIPT_FILES = [
    ("Main_Text_Final_v3.docx", "Main text (SI-synced, submission-ready)"),
    ("Title_Page_Final_v2.docx", "Title page with affiliations"),
    ("Supplementary_Information_Final.docx", "Full supplementary materials (15T + 4F)"),
]


def md5_file(filepath):
    """Calculate MD5 hash of a file."""
    h = hashlib.md5()
    with open(filepath, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            h.update(chunk)
    return h.hexdigest()


def copy_file(src, dst, description=""):
    """Copy a single file with logging."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        shutil.copy2(src, dst)
        size = os.path.getsize(dst)
        return True, size, ""
    except Exception as e:
        return False, 0, str(e)


def main():
    print("=" * 70)
    print("Meniscus_Exo_Multiomics — GitHub Repository Setup Script")
    print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)
    
    stats = {"copied": 0, "skipped": 0, "errors": 0}
    
    # ---- 1. Copy R Scripts ----
    print("\n[1/5] Copying R scripts...")
    for subdir, file_list in R_SCRIPTS.items():
        target_subdir = REPO_ROOT / subdir
        target_subdir.mkdir(parents=True, exist_ok=True)
        
        for item in file_list:
            src_name = item[0]  # relative path or filename
            
            if src_name.startswith("meniscus_exosome_project"):
                src_path = SRC_R_ROOT / src_name
            else:
                src_path = SRC_R_ROOT / src_name
            
            dst_path = target_subdir / Path(src_name).name
            
            if not src_path.exists():
                print(f"  ⚠️  MISSING: {src_path}")
                stats["errors"] += 1
                continue
            
            ok, size, err = copy_file(src_path, dst_path)
            if ok:
                print(f"  ✅ {dst_path.relative_to(REPO_ROOT)} ({size:,} bytes)")
                stats["copied"] += 1
            else:
                print(f"  ❌ {src_path.name}: {err}")
                stats["errors"] += 1
    
    # ---- 2. Copy Python Scripts ----
    print("\n[2/5] Copying Python scripts...")
    for subdir, file_list in PY_SCRIPTS.items():
        target_subdir = REPO_ROOT / subdir
        target_subdir.mkdir(parents=True, exist_ok=True)
        
        for item in file_list:
            src_name = item[0]
            desc = item[1] or ""
            src_path = SRC_PY_ROOT / src_name
            dst_path = target_subdir / src_name
            
            if not src_path.exists():
                print(f"  ⚠️  MISSING: {src_path}")
                stats["errors"] += 1
                continue
            
            ok, size, err = copy_file(src_path, dst_path)
            if ok:
                print(f"  ✅ {dst_path.relative_to(REPO_ROOT)} ({size:,} bytes) | {desc}")
                stats["copied"] += 1
            else:
                print(f"  ❌ {src_name}: {err}")
                stats["errors"] += 1
    
    # ---- 3. Copy CSV Tables ----
    print("\n[3/5] Copying data tables (CSV)...")
    csv_dst_dir = REPO_ROOT / "data" / "tables"
    csv_dst_dir.mkdir(parents=True, exist_ok=True)
    
    # Priority: 1ALL_RESULTS_FINAL first, then fallback to Supplementary_Tables
    csv_src_dir = SRC_DATA_DIR / "results" / "tables" if SRC_DATA_DIR.exists() else SRC_SUPP_TABLES
    
    for csv_name in CSV_FILES:
        src_path = csv_src_dir / csv_name
        
        # Fallback to Supplementary_Tables if not found
        if not src_path.exists():
            alt_src = SRC_SUPP_TABLES / csv_name
            if alt_src.exists():
                src_path = alt_src
        
        if not src_path.exists():
            print(f"  ⚠️  SKIP: {csv_name} (not found)")
            stats["skipped"] += 1
            continue
        
        dst_path = csv_dst_dir / csv_name
        ok, size, err = copy_file(src_path, dst_path)
        if ok:
            print(f"  ✅ data/tables/{csv_name} ({size:,} bytes)")
            stats["copied"] += 1
    
    # ---- 4. Copy Figures ----
    print("\n[4/5] Copying figures...")
    fig_dst_dir = REPO_ROOT / "figures"
    fig_dst_dir.mkdir(parents=True, exist_ok=True)
    
    for fig_name in FIGURE_FILES:
        src_path = SRC_FIGURES / fig_name
        if not src_path.exists():
            print(f"  ⚠️  SKIP: {fig_name} (not found)")
            stats["skipped"] += 1
            continue
        dst_path = fig_dst_dir / fig_name
        ok, size, err = copy_file(src_path, dst_path)
        if ok:
            print(f"  ✅ figures/{fig_name} ({size:,} bytes)")
            stats["copied"] += 1
    
    # ---- 5. Copy Manuscript Files ----
    print("\n[5/5] Copying manuscript files...")
    man_dst_dir = REPO_ROOT / "manuscript"
    man_dst_dir.mkdir(parents=True, exist_ok=True)
    
    for man_item in MANUSCRIPT_FILES:
        man_name = man_item[0]
        desc = man_item[1]
        src_path = SRC_MANUSCRIPT / man_name
        if not src_path.exists():
            print(f"  ⚠️  SKIP: {man_name} (not found)")
            stats["skipped"] += 1
            continue
        dst_path = man_dst_dir / man_name
        ok, size, err = copy_file(src_path, dst_path)
        if ok:
            print(f"  ✅ manuscript/{man_name} ({size:,} bytes) | {desc}")
            stats["copied"] += 1
    
    # ---- Summary ----
    print("\n" + "=" * 70)
    print(f"DONE! Copied: {stats['copied']} | Skipped: {stats['skipped']} | Errors: {stats['errors']}")
    print("=" * 70)
    
    # Generate BOM (Bill of Materials)
    bom_path = REPO_ROOT / "DATA_BOM.md"
    generate_bom(bom_path)
    print(f"\n📋 Data BOM saved to: DATA_BOM.md")


def generate_bom(bom_path):
    """Generate a Bill of Materials markdown file."""
    lines = []
    lines.append("# Data & File Bill of Materials\n")
    lines.append(f"*Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}*\n")
    
    def scan_dir(rel_dir, exts=None):
        abs_dir = REPO_ROOT / rel_dir
        if not abs_dir.exists():
            return []
        items = []
        for f in sorted(abs_dir.iterdir()):
            if f.is_file():
                if exts is None or any(f.name.endswith(e) for e in exts):
                    items.append((f.name, f.stat().st_size))
        return items
    
    # R scripts
    lines.append("## R Analysis Scripts (`R/`)\n")
    lines.append("| File | Size |")
    lines.append("|------|------|")
    total_r_size = 0
    for name, size in scan_dir("R/01_scRNA_integration", [".R"]) + \
                     scan_dir("R/02_exosome_cargo", [".R"]) + \
                     scan_dir("R/03_molecular_docking", [".R"]) + \
                     scan_dir("R/04_validation", [".R"]) + \
                     scan_dir("R/auxiliary", [".R"]):
        lines.append(f"| `R/**/{name}` | {size:,} bytes |")
        total_r_size += size
    lines.append(f"| **Total** | **{total_r_size:,} bytes** |\n")
    
    # Python scripts
    lines.append("## Python Scripts (`python/`)\n")
    lines.append("| File | Size |")
    lines.append("|------|------|")
    total_py_size = 0
    for name, size in scan_dir("python/figures", [".py"]) + \
                     scan_dir("python/analysis", [".py"]) + \
                     scan_dir("python/submission_tools", [".py"]):
        lines.append(f"| `python/**/{name}` | {size:,} bytes |")
        total_py_size += size
    lines.append(f"| **Total** | **{total_py_size:,} bytes** |\n")
    
    # Data tables
    lines.append("## Processed Data Tables (`data/tables/`)\n")
    lines.append("| File | Rows | Size | Description |")
    lines.append("|------|------|------|-------------|")
    import csv
    total_csv_size = 0
    for name, size in scan_dir("data/tables", [".csv"]):
        # Count rows
        try:
            with open(str(REPO_ROOT / "data/tables" / name), 'r', encoding='utf-8') as f:
                reader = csv.reader(f)
                rows = sum(1 for _ in reader) - 1  # minus header
        except:
            rows = "?"
        lines.append(f"| `{name}` | {rows} | {size:,} B | Processed output |")
        total_csv_size += size
    lines.append(f"| **Total** | — | **{total_csv_size:,} B** | — |\n")
    
    # Figures
    lines.append("## Publication Figures (`figures/`)\n")
    lines.append("| File | Size | Format |")
    lines.append("|------|------|--------|")
    total_fig_size = 0
    for name, size in scan_dir("figures"):
        fmt = "Vector PDF" if name.endswith(".pdf") else "Raster PNG (300 DPI)"
        lines.append(f"| `{name}` | {size:,} B | {fmt} |")
        total_fig_size += size
    lines.append(f"| **Total** | **{total_fig_size:,} B** | — |\n")
    
    # Manuscript
    lines.append("## Submission Documents (`manuscript/`)\n")
    lines.append("| Document | Status |")
    lines.append("|----------|--------|")
    for name, size in scan_dir("manuscript"):
        status_map = {
            "Main_Text_Final_v3.docx": "✅ SI-synced (final)",
            "Title_Page_Final_v2.docx": "✅ Canonical author source",
            "Supplementary_Information_Final.docx": "✅ Full version (15T+4F)"
        }
        st = status_map.get(name, "")
        lines.append(f"| `{name}` | {st} |")
    lines.append("")
    
    # Grand total
    grand_total = total_r_size + total_py_size + total_csv_size + total_fig_size
    lines.append("---\n")
    lines.append(f"### Repository Total (excluding .rds): **{grand_total:,} bytes ({grand_total/1024/1024:.1f} MB)**\n")
    
    with open(bom_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))


if __name__ == "__main__":
    main()
