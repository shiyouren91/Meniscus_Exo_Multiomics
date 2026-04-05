# ============================================================
# Step 1: 退变核心轨迹刻画
#
# 核心逻辑: 以正常半月板软骨细胞(FC-ECM)为起点,
#          追踪其向FC-Degenerated / Hypertrophic_Chondrocyte
#          演变的拟时序轨迹, 提取退变驱动基因
#
# 输入: data/processed/meniscus_annotated.rds
# 输出: results/figures/02_trajectory_degeneration/
#       results/tables/trajectory_*.csv
#
# 运行方式: Rscript step1_trajectory_degeneration.R
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(yaml)
  library(tidyverse)
  library(patchwork)
  library(ComplexHeatmap)
  library(ggsci)
  library(grid)
  library(ggrepel)
  library(RColorBrewer)
})

config  <- yaml::read_yaml("config/params.yaml")
data_dir <- path.expand(config$data_dir)
res_dir  <- path.expand(config$results_dir)

fig_dir <- file.path(res_dir, "figures", "02_trajectory_degeneration")
tbl_dir <- file.path(res_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Step 1: 退变核心轨迹刻画\n")
cat("========================================\n\n")

# ============================================================
# 1. 加载数据
# ============================================================
cat("[1/8] 加载Seurat对象...\n")
obj <- readRDS(file.path(data_dir, "processed", "meniscus_annotated.rds"))
cat(paste0("  细胞数: ", ncol(obj), "\n"))
cat(paste0("  条件分组: ", paste(sort(unique(obj$condition)), collapse=", "), "\n"))
cat(paste0("  细胞类型: ", length(unique(obj$cell_type)), " 种\n\n"))

# ============================================================
# 2. 子集: 仅取纤维软骨细胞谱系进行轨迹分析
# ============================================================
cat("[2/8] 子集纤维软骨细胞谱系...\n")

# 定义退变轨迹中的关键细胞类型
trajectory_cell_types <- c(
  "FC-ECM",                          # 正常细胞外基质纤维软骨细胞 (起点)
  "FC-Regulatory",                   # 调控型纤维软骨细胞
  "FC-Metabolic",                    # 代谢型纤维软骨细胞
  "FC-Degenerated",                  # 退变纤维软骨细胞 (终点之一)
  "Hypertrophic_Chondrocyte",        # 肥大软骨细胞 (终点之二)
  "Secretory_Chondrocyte",           # 分泌型软骨细胞
  "Meniscus_Progenitor",             # 半月板祖细胞
  "Stem_Progenitor",                 # 干/祖细胞
  "OuterZone_FC"                     # 外区纤维软骨细胞
)

available_types <- intersect(trajectory_cell_types, unique(obj$cell_type))
missing_types  <- setdiff(trajectory_cell_types, available_types)

if (length(missing_types) > 0) {
  cat(paste0("  警告: 以下细胞类型未找到, 已跳过: ", paste(missing_types, collapse=", "), "\n"))
}

obj_traj <- subset(obj, cell_type %in% available_types)
cat(paste0("  子集后细胞数: ", ncol(obj_traj), "\n"))
cat(paste0("  纳入细胞类型: ", paste(available_types, collapse=", "), "\n\n"))

# ============================================================
# 3. 重新标准化 (使用counts重新处理, 确保Monocle3输入正确)
# ============================================================
cat("[3/8] 数据准备 (SCTransform标准化)...\n")

# 如果没有SCT assay, 使用标准流程
if (!"SCT" %in% names(obj_traj@assays)) {
  obj_traj <- NormalizeData(obj_traj, verbose = FALSE)
  obj_traj <- FindVariableFeatures(obj_traj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  obj_traj <- ScaleData(obj_traj, verbose = FALSE)
  obj_traj <- RunPCA(obj_traj, npcs = 50, verbose = FALSE)
} else {
  DefaultAssay(obj_traj) <- "SCT"
}

obj_traj <- RunUMAP(obj_traj, reduction = "pca", dims = 1:30, verbose = FALSE)

cat("  UMAP计算完成\n\n")

# ============================================================
# 4. Monocle3 拟时序分析
# ============================================================
cat("[4/8] Monocle3 拟时序分析...\n")

# Seurat → cell_data_set 转换
expr_mat <- GetAssayData(obj_traj, slot = "counts")
cell_meta <- obj_traj@meta.data
gene_meta <- data.frame(gene_short_name = rownames(expr_mat),
                        row.names = rownames(expr_mat))

cds <- new_cell_data_set(expr_mat,
                         cell_metadata = cell_meta,
                         gene_metadata = gene_meta)

# 预处理 + 导入Seurat UMAP
cds <- preprocess_cds(cds, num_dim = 30)
reducedDims(cds)[["UMAP"]] <- Embeddings(obj_traj, "umap")

# 聚类 (使用Seurat已有的细胞类型标签)
cds <- cluster_cells(cds, reduction_method = "UMAP", resolution = 1e-3)  # 极低resolution合并为大cluster

# 学习轨迹图
cds <- learn_graph(cds, use_partition = TRUE)

# ============================================================
# 5. 设定根节点: FC-ECM细胞 (正常半月板软骨细胞)
# ============================================================
cat("[5/8] 设定根节点 (FC-ECM = 退变起点)...\n")

# FC-ECM的marker基因
ecm_markers <- c("COL1A1", "COL2A1", "ACAN", "COMP", "PRG4", "LUM", "DCN")
ecm_markers <- ecm_markers[ecm_markers %in% rownames(cds)]

# 计算每个细胞的FC-ECM特征评分
if (length(ecm_markers) > 0) {
  ecm_score <- colMeans(exprs(cds)[ecm_markers, , drop = FALSE])
  root_cells <- names(sort(ecm_score, decreasing = TRUE))[1:5]
  cds <- order_cells(cds, root_cells = root_cells)
  cat(paste0("  根节点细胞数: ", length(root_cells), "\n"))
} else {
  # 如果marker不存在, 直接用cell_type选择FC-ECM细胞
  ecm_cells <- colnames(obj_traj)[obj_traj$cell_type == "FC-ECM"]
  if (length(ecm_cells) > 0) {
    cds <- order_cells(cds, root_cells = ecm_cells[1:min(5, length(ecm_cells))])
  } else {
    cds <- order_cells(cds)
    cat("  警告: 未找到FC-ECM细胞, 使用默认根节点\n")
  }
}

# 保存pseudotime到Seurat对象
obj_traj$pseudotime <- pseudotime(cds)
obj_traj$principal_graph <- as.character(cds@principal_graph_aux$principal_graph[[1]]$df$cell)

cat("  拟时序计算完成\n\n")

# ============================================================
# 6. 提取退变驱动基因
# ============================================================
cat("[6/8] 提取退变驱动基因...\n")

# 6.1 Monocle3 graph_test: 沿轨迹显著变化的基因
trajectory_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 1)
trajectory_genes <- trajectory_genes %>%
  filter(q_value < 0.05) %>%
  arrange(q_value)

cat(paste0("  沿轨迹显著变化基因 (q<0.05): ", nrow(trajectory_genes), "\n"))

# 6.2 按pseudotime分bin, 做时间序列差异分析
# 将pseudotime分为4个阶段
pt <- obj_traj$pseudotime
pt[is.na(pt)] <- 0
pt_breaks <- quantile(pt, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
pt_breaks <- unique(pt_breaks)
obj_traj$pseudotime_bin <- cut(pt, breaks = pt_breaks, 
                                labels = c("Early", "Mid-Early", "Mid-Late", "Late"),
                                include.lowest = TRUE)

# Late vs Early 差异基因 (退变终末期 vs 正常早期)
DefaultAssay(obj_traj) <- "RNA"
Idents(obj_traj) <- "pseudotime_bin"
degs_late_vs_early <- tryCatch({
  FindMarkers(obj_traj, ident.1 = "Late", ident.2 = "Early",
              min.pct = 0.1, logfc.threshold = 0.25) %>%
    rownames_to_column(var = "gene") %>%
    filter(p_val_adj < 0.05)
}, error = function(e) {
  cat(paste0("  警告: 差异分析失败: ", e$message, "\n"))
  data.frame(gene = character(), stringsAsFactors = FALSE)
})

if (nrow(degs_late_vs_early) > 0) {
  cat(paste0("  Late vs Early 显著差异基因: ", nrow(degs_late_vs_early), "\n"))
}

# 6.3 合并为"退变驱动基因"
driver_genes <- trajectory_genes %>%
  rownames_to_column(var = "gene") %>%
  left_join(degs_late_vs_early %>% select(gene, avg_log2FoldChange = avg_log2FoldChange),
            by = "gene")

# 标注方向: Late期上调 = 退变促进基因
if ("avg_log2FoldChange" %in% colnames(driver_genes)) {
  driver_genes <- driver_genes %>%
    mutate(direction = case_when(
      avg_log2FoldChange > 0.5 ~ "Degeneration_UP",
      avg_log2FoldChange < -0.5 ~ "Degeneration_DOWN",
      TRUE ~ "Neutral"
    ))
}

write.csv(driver_genes,
          file.path(tbl_dir, "trajectory_degeneration_driver_genes.csv"),
          row.names = FALSE)

cat(paste0("  退变驱动基因已保存\n\n"))

# ============================================================
# 7. 绘制论文级图表
# ============================================================
cat("[7/8] 绘制论文级图表...\n")

# --- 7.1 核心图: 拟时序轨迹 + 细胞类型着色 ---
pdf(file.path(fig_dir, "Fig_trajectory_overview.pdf"), width = 18, height = 7)
p1 <- plot_cells(cds, color_cells_by = "pseudotime",
                 cell_size = 0.4, label_branch_points = TRUE, label_roots = TRUE,
                 show_trajectory_graph = TRUE) +
  scale_color_viridis_c(name = "Pseudotime", option = "plasma") +
  ggtitle("Degeneration Trajectory") +
  theme(legend.position = "right")

p2 <- plot_cells(cds, color_cells_by = "cell_type",
                 cell_size = 0.4, label_groups_by_cluster = TRUE,
                 group_label_size = 3.5, show_trajectory_graph = TRUE) +
  ggtitle("Trajectory by Cell Type")

print(p1 | p2)
dev.off()

# --- 7.2 UMAP 按pseudotime_bin着色 ---
pdf(file.path(fig_dir, "Fig_trajectory_pseudotime_bins.pdf"), width = 10, height = 8)
print(DimPlot(obj_traj, group.by = "pseudotime_bin", 
              cols = c("Early" = "#4DBBD5", "Mid-Early" = "#91D1C2",
                       "Mid-Late" = "#FDB863", "Late" = "#E64B35"),
              pt.size = 0.4, order = c("Early", "Mid-Early", "Mid-Late", "Late")) +
      ggtitle("Pseudotime Stages: Early → Late (Degeneration)"))
dev.off()

# --- 7.3 UMAP 按condition着色 (确认轨迹与退变状态的对应关系) ---
pdf(file.path(fig_dir, "Fig_trajectory_condition.pdf"), width = 14, height = 7)
print(DimPlot(obj_traj, group.by = "condition", split.by = "condition",
              cols = c("Normal" = "#4DBBD5", "Degenerated" = "#E64B35"),
              pt.size = 0.3) +
      ggtitle("Normal vs Degenerated on Trajectory UMAP"))
dev.off()

# --- 7.4 关键基因沿拟时序的表达动态 (论文核心图) ---
key_genes <- c(
  # ECM合成基因 (早期高, 晚期低)
  "COL1A1", "COL2A1", "ACAN", "PRG4", "COMP", "DCN", "LUM",
  # 退变相关基因 (晚期高)
  "MMP1", "MMP3", "MMP13", "ADAMTS4", "ADAMTS5",
  # 炎症基因 (晚期高)
  "IL1B", "IL6", "TNF", "PTGS2", "CXCL8",
  # 肥大标志 (晚期高)
  "COL10A1", "RUNX2", "ALPL", "IHH",
  # 修复/再生相关
  "SOX9", "TGFB1", "TGFB3", "BMP2", "IGF1",
  # 凋亡/衰老
  "CDKN1A", "CDKN2A", "BAX", "BCL2"
)
key_genes <- key_genes[key_genes %in% rownames(cds)]

if (length(key_genes) > 0) {
  pdf(file.path(fig_dir, "Fig_trajectory_gene_dynamic.pdf"), width = 24, height = 20)
  plot_genes_in_pseudotime(
    cds[key_genes, ],
    color_cells_by = "cell_type",
    min_expr = 0.5,
    ncol = 6,
    alpha = 0.3
  )
  dev.off()
}

# --- 7.5 退变驱动基因热图 (按pseudotime阶段聚类) ---
if (nrow(driver_genes) > 0) {
  top_driver_genes <- head(driver_genes$gene, 50)
  top_driver_genes <- top_driver_genes[top_driver_genes %in% rownames(obj_traj)]
  
  if (length(top_driver_genes) > 0) {
    # 按pseudotime分bin的平均表达
    bins_order <- c("Early", "Mid-Early", "Mid-Late", "Late")
    bins_order <- intersect(bins_order, levels(obj_traj$pseudotime_bin))
    
    bin_expr <- AverageExpression(obj_traj, features = top_driver_genes,
                                   group.by = "pseudotime_bin")$RNA
    bin_expr <- bin_expr[, bins_order, drop = FALSE]
    
    # Z-score标准化
    bin_expr_z <- t(scale(t(as.matrix(bin_expr))))
    bin_expr_z[is.na(bin_expr_z)] <- 0
    
    # 注释条: 基因方向
    direction_anno <- driver_genes %>%
      filter(gene %in% top_driver_genes) %>%
      select(gene, direction)
    direction_anno <- direction_anno$direction[match(top_driver_genes, direction_anno$gene)]
    direction_anno[is.na(direction_anno)] <- "Neutral"
    
    anno_row <- rowAnnotation(
      Direction = direction_anno,
      col = list(Direction = c("Degeneration_UP" = "#E64B35",
                                "Degeneration_DOWN" = "#4DBBD5",
                                "Neutral" = "grey70"))
    )
    
    pdf(file.path(fig_dir, "Fig_trajectory_driver_genes_heatmap.pdf"), width = 10, height = 14)
    Heatmap(
      bin_expr_z,
      name = "Z-score",
      col = colorRamp2(c(-2, 0, 2), c("#4DBBD5", "white", "#E64B35")),
      cluster_rows = TRUE,
      show_row_names = TRUE,
      row_names_gp = gpar(fontsize = 7),
      column_names_rot = 45,
      column_title = "Early → Late (Degeneration)",
      row_title = "Top 50 Driver Genes",
      right_annotation = anno_row,
      heatmap_legend_param = list(title_gp = gpar(fontsize = 10))
    )
    dev.off()
  }
}

# --- 7.6 各细胞类型沿pseudotime的分布 (箱线图) ---
pdf(file.path(fig_dir, "Fig_pseudotime_by_celltype.pdf"), width = 14, height = 6)
cell_order <- c("FC-ECM", "FC-Regulatory", "Meniscus_Progenitor", "Stem_Progenitor",
                "OuterZone_FC", "FC-Metabolic", "Secretory_Chondrocyte",
                "FC-Degenerated", "Hypertrophic_Chondrocyte")
cell_order <- intersect(cell_order, available_types)

p_box <- ggplot(obj_traj@meta.data, aes(x = factor(cell_type, levels = cell_order), 
                                          y = pseudotime, fill = cell_type)) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.3) +
  scale_fill_npg() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none") +
  labs(title = "Pseudotime Distribution by Cell Type",
       subtitle = "FC-ECM → FC-Degenerated / Hypertrophic_Chondrocyte",
       x = "", y = "Pseudotime")
print(p_box)
dev.off()

# --- 7.7 条件(正常vs退变)沿pseudotime的分布 ---
pdf(file.path(fig_dir, "Fig_pseudotime_by_condition.pdf"), width = 10, height = 5)
p_cond <- ggplot(obj_traj@meta.data, aes(x = condition, y = pseudotime, fill = condition)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.size = 0.3) +
  scale_fill_manual(values = c("Normal" = "#4DBBD5", "Degenerated" = "#E64B35")) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Normal vs Degenerated: Pseudotime Distribution",
       x = "", y = "Pseudotime")

# 统计检验
wilcox_test <- wilcox.test(
  obj_traj$pseudotime[obj_traj$condition == "Degenerated"],
  obj_traj$pseudotime[obj_traj$condition == "Normal"]
)
p_val_label <- ifelse(wilcox_test$p.value < 0.001, 
                       paste0("p < 0.001\n(Wilcoxon)"),
                       paste0("p = ", format(wilcox_test$p.value, digits = 3), "\n(Wilcoxon)"))
p_cond <- p_cond + annotate("text", x = 1.5, y = max(obj_traj$pseudotime, na.rm = TRUE) * 0.95,
                              label = p_val_label, size = 5, fontface = "bold")
print(p_cond)
dev.off()

cat(paste0("  Wilcoxon test: p = ", format(wilcox_test$p.value, digits = 4), "\n\n"))

# --- 7.8 退变路径Sankey图: 细胞状态转换流 ---
sankey_data <- obj_traj@meta.data %>%
  group_by(pseudotime_bin, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(pseudotime_bin) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

pdf(file.path(fig_dir, "Fig_trajectory_sankey.pdf"), width = 14, height = 8)
p_sankey <- ggplot(sankey_data,
                   aes(axis1 = pseudotime_bin, axis2 = cell_type, y = prop)) +
  geom_alluvium(aes(fill = cell_type), alpha = 0.7, width = 1/12) +
  geom_stratum(width = 1/8, fill = "grey90", color = "grey30", size = 0.3) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
  scale_fill_npg() +
  theme_minimal() +
  labs(title = "Cell State Transition Along Degeneration Trajectory",
       subtitle = "FC-ECM (Normal) → FC-Degenerated / Hypertrophic (Late)",
       y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        legend.position = "none",
        axis.title.x = element_blank())
print(p_sankey)
dev.off()

# --- 7.9 轨迹分化分支鉴定 ---
# 提取分支信息
if (!is.null(cds@principal_graph_aux)) {
  # 获取每个细胞所在的branch
  # Monocle3的默认branch: 每个cell的closest_dp
  closest_dp <- cds@principal_graph_aux$dp_mats[[1]]
  
  # 如果能获取branch分配
  tryCatch({
    pdf(file.path(fig_dir, "Fig_trajectory_branches.pdf"), width = 14, height = 7)
    p_branch <- plot_cells(cds, color_cells_by = "cell_type",
                           cell_size = 0.3,
                           label_branch_points = TRUE,
                           label_roots = TRUE,
                           show_trajectory_graph = TRUE) +
      ggtitle("Degeneration Trajectory Branches") +
      theme(legend.position = "right",
            legend.text = element_text(size = 8))
    print(p_branch)
    dev.off()
  }, error = function(e) {
    cat(paste0("  Branch plot skipped: ", e$message, "\n"))
  })
}

# ============================================================
# 8. 保存结果
# ============================================================
cat("[8/8] 保存结果...\n")

# 更新Seurat对象
saveRDS(obj_traj, file.path(data_dir, "processed", "meniscus_trajectory_degeneration.rds"))

# 保存CDS对象 (Monocle3)
saveRDS(cds, file.path(data_dir, "processed", "cds_degeneration.rds"))

# 保存驱动基因统计
gene_stats <- data.frame(
  Total_cells = ncol(obj_traj),
  Cell_types_used = length(available_types),
  Trajectory_genes_sig = nrow(driver_genes),
  Degeneration_UP_genes = sum(driver_genes$direction == "Degeneration_UP", na.rm = TRUE),
  Degeneration_DOWN_genes = sum(driver_genes$direction == "Degeneration_DOWN", na.rm = TRUE),
  Wilcoxon_pvalue = wilcox_test$p.value
)
write.csv(gene_stats, file.path(tbl_dir, "trajectory_analysis_summary.csv"), row.names = FALSE)

# 打印核心结果摘要
cat("\n========================================\n")
cat("分析完成! 核心结果摘要:\n")
cat("========================================\n")
cat(paste0("  纳入细胞: ", ncol(obj_traj), " (来自", length(available_types), "种细胞类型)\n"))
cat(paste0("  轨迹显著变化基因: ", nrow(driver_genes), "\n"))
cat(paste0("  退变上调基因: ", sum(driver_genes$direction == "Degeneration_UP", na.rm = TRUE), "\n"))
cat(paste0("  退变下调基因: ", sum(driver_genes$direction == "Degeneration_DOWN", na.rm = TRUE), "\n"))
cat(paste0("  Normal vs Deg pseudotime差异: p = ", format(wilcox_test$p.value, digits = 4), "\n"))
cat("\n输出文件:\n")
cat(paste0("  图表: ", fig_dir, "\n"))
cat(paste0("  表格: ", tbl_dir, "/trajectory_*.csv\n"))
cat(paste0("  Seurat对象: data/processed/meniscus_trajectory_degeneration.rds\n"))
cat(paste0("  CDS对象: data/processed/cds_degeneration.rds\n"))
cat("\n下一步: Rscript step2_cellchat_differential.R\n")
