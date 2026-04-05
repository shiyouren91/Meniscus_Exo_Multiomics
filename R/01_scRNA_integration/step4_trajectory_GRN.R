# ============================================================
# Step 4: 拟时序分析与转录调控网络
# 
# 输入: 注释完成的Seurat对象
# 输出: 拟时序轨迹 + SCENIC调控网络
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(slingshot)
  library(yaml)
  library(tidyverse)
  library(patchwork)
  library(ComplexHeatmap)
  library(ggsci)
})

config  <- yaml::read_yaml("config/params.yaml")
data_dir <- path.expand(config$data_dir)
res_dir  <- path.expand(config$results_dir)

dir.create(file.path(res_dir, "figures", "02_trajectory"), recursive = TRUE, showWarnings = FALSE)

cat("=== Step 4: Trajectory & GRN Analysis ===\n")

merged <- readRDS(file.path(data_dir, "processed", "meniscus_annotated.rds"))

# ============================================================
# 4.1 Monocle3 拟时序分析
# ============================================================

cat("  Running Monocle3 pseudotime analysis...\n")

# Seurat -> Monocle3 转换
expression_matrix <- GetAssayData(merged, slot = "counts")
cell_metadata <- merged@meta.data
gene_metadata <- data.frame(
  gene_short_name = rownames(merged),
  row.names = rownames(merged)
)

cds <- new_cell_data_set(
  expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

# 使用Seurat的UMAP坐标
cds <- preprocess_cds(cds, num_dim = 30)

# 导入Seurat UMAP
umap_coords <- Embeddings(merged, "umap")
reducedDims(cds)[["UMAP"]] <- umap_coords

# 聚类（使用Seurat的聚类结果）
cds <- cluster_cells(cds, reduction_method = "UMAP")

# 学习轨迹图
cds <- learn_graph(cds, use_partition = TRUE)

# 选择根节点（以祖细胞群/EC为起点）
# 自动选择：找到表达祖细胞marker最高的节点
progenitor_markers <- c("MCAM", "PDGFRB", "ENG", "CD34")
progenitor_markers <- progenitor_markers[progenitor_markers %in% rownames(cds)]

if (length(progenitor_markers) > 0) {
  # 计算每个细胞的祖细胞marker平均表达
  prog_expr <- colMeans(exprs(cds)[progenitor_markers, , drop = FALSE])
  
  # 选择表达最高的细胞作为根
  root_cells <- names(sort(prog_expr, decreasing = TRUE))[1:10]
  
  cds <- order_cells(cds, root_cells = root_cells)
} else {
  cat("  Warning: No progenitor markers found. Manual root selection needed.\n")
  # 使用partition中最早的点
  cds <- order_cells(cds)
}

# 拟时序可视化
p_pseudo1 <- plot_cells(cds, color_cells_by = "pseudotime",
                        cell_size = 0.5, label_groups_by_cluster = FALSE,
                        label_branch_points = TRUE, label_roots = TRUE) +
  ggtitle("Pseudotime Trajectory")

p_pseudo2 <- plot_cells(cds, color_cells_by = "cell_type",
                        cell_size = 0.5, label_groups_by_cluster = TRUE,
                        group_label_size = 3) +
  ggtitle("Trajectory colored by Cell Type")

pdf(file.path(res_dir, "figures", "02_trajectory", "monocle3_pseudotime.pdf"),
    width = 16, height = 7)
p_pseudo1 + p_pseudo2
dev.off()

# 保存pseudotime到Seurat对象
merged$pseudotime <- pseudotime(cds)

# ============================================================
# 4.2 沿拟时序的基因表达动态
# ============================================================

cat("  Analyzing gene dynamics along pseudotime...\n")

# 找到沿轨迹变化的基因
trajectory_genes <- graph_test(cds, neighbor_graph = "principal_graph",
                                cores = 4)

trajectory_genes <- trajectory_genes %>%
  filter(q_value < 0.05) %>%
  arrange(q_value)

write.csv(trajectory_genes,
          file.path(res_dir, "tables", "trajectory_dynamic_genes.csv"),
          row.names = TRUE)

# Top动态基因可视化
top_trajectory_genes <- head(rownames(trajectory_genes), 20)

pdf(file.path(res_dir, "figures", "02_trajectory", "trajectory_gene_expression.pdf"),
    width = 20, height = 20)
plot_cells(cds, genes = top_trajectory_genes, cell_size = 0.3,
           show_trajectory_graph = FALSE) &
  scale_color_viridis_c()
dev.off()

# 治疗相关基因沿拟时序的表达
therapeutic_genes <- c(
  # ECM 合成
  "COL1A1", "COL2A1", "ACAN", "SOX9", "PRG4",
  # 基质降解
  "MMP1", "MMP3", "MMP13", "ADAMTS4", "ADAMTS5",
  # 炎症
  "IL1B", "IL6", "TNF", "PTGS2", "CCL2",
  # 生长因子/信号
  "TGFB1", "TGFB3", "FGF2", "IGF1", "BMP2",
  # 凋亡/衰老
  "CDKN1A", "CDKN2A", "TP53", "BCL2", "BAX"
)
therapeutic_genes <- therapeutic_genes[therapeutic_genes %in% rownames(cds)]

pdf(file.path(res_dir, "figures", "02_trajectory", "therapeutic_genes_pseudotime.pdf"),
    width = 24, height = 20)
plot_genes_in_pseudotime(
  cds[therapeutic_genes, ],
  color_cells_by = "cell_type",
  min_expr = 0.5,
  ncol = 4
)
dev.off()

# ============================================================
# 4.3 Slingshot 轨迹分析（补充验证）
# ============================================================

cat("  Running Slingshot trajectory analysis...\n")

# 提取UMAP坐标和聚类信息
umap_emb <- Embeddings(merged, "umap")
cluster_labels <- merged$cell_type

# Slingshot
sds <- slingshot(umap_emb, clusterLabels = cluster_labels)

# 保存slingshot结果
merged$slingshot_pseudotime <- slingPseudotime(sds)[, 1]  # 第一条lineage

# Slingshot 可视化
pdf(file.path(res_dir, "figures", "02_trajectory", "slingshot_trajectory.pdf"),
    width = 10, height = 8)
plot(umap_emb, col = alpha(colors()[sample(1:600, length(unique(cluster_labels)))][as.factor(cluster_labels)], 0.5),
     pch = 16, cex = 0.3, xlab = "UMAP1", ylab = "UMAP2",
     main = "Slingshot Trajectories")
lines(SlingshotDataSet(sds), lwd = 2, type = "lineages", col = "black")
dev.off()

# ============================================================
# 4.4 CytoTRACE 分化潜能分析
# ============================================================

cat("  Running CytoTRACE analysis...\n")

# CytoTRACE 评估细胞分化状态
# 基于基因数量的简化CytoTRACE计算
cytotrace_score <- colSums(GetAssayData(merged, slot = "counts") > 0)
merged$cytotrace_score <- cytotrace_score / max(cytotrace_score)

p_cyto <- FeaturePlot(merged, features = "cytotrace_score", 
                       cols = c("blue", "yellow", "red")) +
  ggtitle("CytoTRACE Score (Differentiation Potential)")

pdf(file.path(res_dir, "figures", "02_trajectory", "cytotrace_score.pdf"),
    width = 10, height = 8)
print(p_cyto)
dev.off()

# 按细胞类型比较分化潜能
p_cyto_box <- VlnPlot(merged, features = "cytotrace_score",
                       group.by = "cell_type", pt.size = 0) +
  RotatedAxis() +
  ggtitle("CytoTRACE by Cell Type")

pdf(file.path(res_dir, "figures", "02_trajectory", "cytotrace_by_celltype.pdf"),
    width = 12, height = 5)
print(p_cyto_box)
dev.off()

# ============================================================
# 4.5 SCENIC 转录因子调控网络（简化版）
# ============================================================

cat("  Preparing SCENIC analysis...\n")
cat("  Note: Full SCENIC requires pySCENIC or R SCENIC pipeline.\n")
cat("  Generating input files for SCENIC...\n")

# 导出SCENIC所需的输入文件
scenic_dir <- file.path(data_dir, "processed", "scenic_input")
dir.create(scenic_dir, recursive = TRUE, showWarnings = FALSE)

# 表达矩阵（loom格式需要另外转换，这里导出为csv）
# 取子集以减少计算量（随机采样5000细胞）
set.seed(42)
cells_subset <- sample(colnames(merged), min(5000, ncol(merged)))
expr_matrix <- GetAssayData(merged, slot = "counts")[, cells_subset]

# 保存为Matrix Market格式
Matrix::writeMM(expr_matrix, file.path(scenic_dir, "expression_matrix.mtx"))
write.csv(rownames(expr_matrix), file.path(scenic_dir, "genes.csv"), row.names = FALSE)
write.csv(colnames(expr_matrix), file.path(scenic_dir, "cells.csv"), row.names = FALSE)
write.csv(merged@meta.data[cells_subset, c("cell_type", "dataset")],
          file.path(scenic_dir, "cell_metadata.csv"))

# 生成pySCENIC运行脚本
cat('#!/bin/bash
# ============================================================
# pySCENIC Pipeline
# 运行前确保已安装 pyscenic 和下载参考数据库
# ============================================================

# 参考数据库下载 (约2GB)
# wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
# wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10-nr.hgnc-m0.001-o0.0.tbl

SCENIC_DIR="scenic_input"
DB_DIR="scenic_databases"  # 修改为你的数据库路径

# Step 1: GRN inference
pyscenic grn \\
  ${SCENIC_DIR}/expression_matrix.loom \\
  ${DB_DIR}/allTFs_hg38.txt \\
  -o ${SCENIC_DIR}/adjacencies.csv \\
  --num_workers 8

# Step 2: Regulon prediction (cisTarget)
pyscenic ctx \\
  ${SCENIC_DIR}/adjacencies.csv \\
  ${DB_DIR}/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \\
  --annotations_fname ${DB_DIR}/motifs-v10-nr.hgnc-m0.001-o0.0.tbl \\
  --expression_mtx_fname ${SCENIC_DIR}/expression_matrix.loom \\
  --output ${SCENIC_DIR}/regulons.csv \\
  --num_workers 8

# Step 3: AUCell scoring
pyscenic aucell \\
  ${SCENIC_DIR}/expression_matrix.loom \\
  ${SCENIC_DIR}/regulons.csv \\
  --output ${SCENIC_DIR}/auc_matrix.csv \\
  --num_workers 8

echo "SCENIC analysis complete!"
', file = file.path(scenic_dir, "run_pyscenic.sh"))

# ============================================================
# 4.6 关键转录因子与外泌体cargo的关联分析
# ============================================================

cat("  Analyzing TF-exosome cargo connections...\n")

# 已知与半月板/软骨修复相关的关键转录因子
key_tfs <- c("SOX9", "SOX5", "SOX6",      # 软骨分化
             "RUNX2", "RUNX1",              # 骨/软骨分化
             "SMAD2", "SMAD3", "SMAD4",     # TGF-β信号
             "NFKB1", "RELA",               # NF-κB炎症
             "STAT3", "STAT1",              # JAK/STAT
             "HIF1A",                        # 缺氧响应
             "FOXO1", "FOXO3",              # 凋亡/自噬
             "MYC", "JUN", "FOS",           # 增殖/AP-1
             "ETS1", "ETS2",                # 血管生成
             "ZEB1", "TWIST1")              # EMT/衰老

key_tfs <- key_tfs[key_tfs %in% rownames(merged)]

# TF在各细胞类型中的表达
tf_expr <- AverageExpression(merged, features = key_tfs, group.by = "cell_type")$RNA

# 热图
pdf(file.path(res_dir, "figures", "02_trajectory", "key_TF_heatmap.pdf"),
    width = 10, height = 8)
Heatmap(
  as.matrix(tf_expr),
  name = "Expression",
  col = circlize::colorRamp2(c(0, 1, 3), c("blue", "white", "red")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 45,
  row_title = "Transcription Factors",
  column_title = "Cell Types"
)
dev.off()

# ============================================================
# 4.7 保存
# ============================================================

saveRDS(merged, file.path(data_dir, "processed", "meniscus_trajectory.rds"))

cat("\n  Trajectory analysis complete.\n")
cat("  SCENIC input files saved to: scenic_input/\n")
cat("  Next step: Rscript scripts/02_exosome_cargo/step1_cargo_database.R\n")
cat("=== Step 4 Complete ===\n")
