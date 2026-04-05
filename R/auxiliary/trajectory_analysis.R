# ============================================================
# 轨迹分析：正常 → 退变的细胞状态转换
# 上传到 AutoDL 后执行:
#   cd /root/autodl-tmp/meniscus_exosome_project
#   Rscript trajectory_analysis.R 2>&1 | tee logs/trajectory.log
# ============================================================

options(conflicts.policy = list(warn = FALSE))
suppressPackageStartupMessages({
  library(Seurat)
  library(slingshot)
  library(ggplot2)
  library(dplyr)
  library(yaml)
  library(ComplexHeatmap)
  library(circlize)
})

config <- yaml::read_yaml("config/params.yaml")
data_dir <- config$data_dir; res_dir <- config$results_dir
dir.create(file.path(res_dir, "figures/02_trajectory"), recursive = TRUE, showWarnings = FALSE)

cat("=== Trajectory Analysis ===\n")
obj <- readRDS(file.path(data_dir, "processed/meniscus_annotated.rds"))
cat(paste0("  Cells: ", ncol(obj), " | Cell types: ", length(unique(obj$cell_type)), "\n"))

# ============================================================
# 1. Slingshot 轨迹推断
# ============================================================
cat("\n--- 1. Slingshot Trajectory ---\n")

# 提取纤维软骨细胞亚群（核心退变轨迹）
fc_types <- c("FC-Regulatory", "FC-ECM", "FC-Metabolic", "FC-Inflammatory",
              "FC-Degenerated", "Secretory_Chondrocyte", "Hypertrophic_Chondrocyte",
              "OuterZone_FC", "Meniscus_Progenitor")
fc_cells <- subset(obj, cell_type %in% fc_types)
cat(paste0("  FC subset: ", ncol(fc_cells), " cells\n"))

# 用UMAP坐标做Slingshot
umap_emb <- Embeddings(fc_cells, "umap")
cl <- fc_cells$cell_type

sds <- slingshot(umap_emb, clusterLabels = cl, start.clus = "FC-Regulatory")

# 提取pseudotime
pt <- slingPseudotime(sds)
n_lineages <- ncol(pt)
cat(paste0("  Lineages found: ", n_lineages, "\n"))

# 把pseudotime存入Seurat对象
for (i in 1:n_lineages) {
  col_name <- paste0("pseudotime_", i)
  fc_cells@meta.data[[col_name]] <- pt[, i]
}

# 主pseudotime（取lineage 1）
fc_cells$pseudotime <- pt[, 1]

# 也存回全部对象
obj@meta.data$pseudotime <- NA
obj@meta.data[colnames(fc_cells), "pseudotime"] <- fc_cells$pseudotime

# ============================================================
# 2. 轨迹可视化
# ============================================================
cat("\n--- 2. Trajectory Visualization ---\n")

# 2a. UMAP + 轨迹曲线
colors_fc <- c("FC-Regulatory"="#E64B35", "FC-ECM"="#F39B7F", "FC-Metabolic"="#E7298A",
               "FC-Inflammatory"="#E6AB02", "FC-Degenerated"="#7570B3",
               "Secretory_Chondrocyte"="#D95F02", "Hypertrophic_Chondrocyte"="#FB8072",
               "OuterZone_FC"="#FDB462", "Meniscus_Progenitor"="#1B9E77")

# 提取Slingshot曲线坐标
curves_data <- list()
for (i in seq_along(slingCurves(sds))) {
  crv <- slingCurves(sds)[[i]]$s
  curves_data[[i]] <- data.frame(x = crv[, 1], y = crv[, 2], lineage = i)
}
curves_df <- do.call(rbind, curves_data)

plot_df <- data.frame(
  UMAP1 = umap_emb[, 1], UMAP2 = umap_emb[, 2],
  cell_type = cl, pseudotime = fc_cells$pseudotime
)

p_traj <- ggplot(plot_df, aes(UMAP1, UMAP2)) +
  geom_point(aes(color = cell_type), size = 0.5, alpha = 0.6) +
  geom_path(data = curves_df, aes(x, y, group = lineage), linewidth = 1.2, color = "black") +
  scale_color_manual(values = colors_fc) +
  theme_classic() +
  labs(title = "Slingshot Trajectory: FC Subtypes", color = "Cell Type") +
  theme(legend.text = element_text(size = 8))

pdf(file.path(res_dir, "figures/02_trajectory/slingshot_trajectory.pdf"), width = 12, height = 8)
print(p_traj)
dev.off()

# 2b. Pseudotime 着色
p_pt <- ggplot(plot_df, aes(UMAP1, UMAP2, color = pseudotime)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c(option = "inferno", na.value = "grey90") +
  theme_classic() +
  labs(title = "Pseudotime (FC-Regulatory → Degeneration)", color = "Pseudotime")

pdf(file.path(res_dir, "figures/02_trajectory/pseudotime_umap.pdf"), width = 10, height = 8)
print(p_pt)
dev.off()

# 2c. 条件沿pseudotime分布
plot_df$condition <- fc_cells$condition
p_cond_pt <- ggplot(plot_df[!is.na(plot_df$pseudotime), ], aes(pseudotime, fill = condition)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("Normal" = "#4DBBD5", "Degenerated" = "#E64B35")) +
  theme_classic() +
  labs(title = "Pseudotime Distribution: Normal vs Degenerated",
       x = "Pseudotime", y = "Density")

pdf(file.path(res_dir, "figures/02_trajectory/pseudotime_condition_density.pdf"), width = 8, height = 5)
print(p_cond_pt)
dev.off()

# ============================================================
# 3. 沿轨迹的关键基因表达动态
# ============================================================
cat("\n--- 3. Gene Dynamics Along Pseudotime ---\n")

# 按pseudotime排序，对关键基因做平滑
key_genes <- c(
  # 正常标志
  "COL1A1", "COL2A1", "ACAN", "SOX9", "PRG4",
  # 退变标志
  "MMP13", "ADAMTS5", "GDF15", "CD274",
  # 炎症
  "IL33", "IL1RL1", "IL1B", "TNF",
  # 修复/外泌体靶向
  "TGFB1", "TGFB3", "FGF2", "MCAM",
  # 信号受体
  "TGFBR1", "TGFBR2", "FGFR1", "BMPR1A"
)
key_genes <- key_genes[key_genes %in% rownames(fc_cells)]

# 获取表达数据
expr_data <- GetAssayData(fc_cells, layer = "data")[key_genes, ]
meta <- fc_cells@meta.data
valid_cells <- !is.na(meta$pseudotime)
expr_valid <- as.matrix(expr_data[, valid_cells])
pt_valid <- meta$pseudotime[valid_cells]

# 按pseudotime排序
ord <- order(pt_valid)
expr_ordered <- expr_valid[, ord]
pt_ordered <- pt_valid[ord]

# 滑动窗口平滑
window_size <- min(100, ncol(expr_ordered) %/% 10)
smoothed <- t(apply(expr_ordered, 1, function(x) {
  stats::filter(x, rep(1/window_size, window_size), sides = 2)
}))
colnames(smoothed) <- 1:ncol(smoothed)

# 去掉NA列
na_cols <- apply(smoothed, 2, function(x) any(is.na(x)))
smoothed_clean <- smoothed[, !na_cols]

# Z-score标准化
smoothed_z <- t(scale(t(smoothed_clean)))
smoothed_z[is.na(smoothed_z)] <- 0

# 热图
n_cols <- ncol(smoothed_z)
col_split <- cut(1:n_cols, breaks = 5, labels = c("Early", "Early-Mid", "Mid", "Mid-Late", "Late"))

pdf(file.path(res_dir, "figures/02_trajectory/gene_dynamics_heatmap.pdf"), width = 14, height = 10)
Heatmap(smoothed_z,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("#4DBBD5", "white", "#E64B35")),
  cluster_rows = TRUE, cluster_columns = FALSE,
  show_column_names = FALSE, show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_title = "Pseudotime (Early → Late Degeneration)",
  row_title = "Key Genes",
  top_annotation = HeatmapAnnotation(
    Stage = col_split,
    col = list(Stage = c("Early" = "#4DBBD5", "Early-Mid" = "#91D1C2",
                         "Mid" = "#F0E442", "Mid-Late" = "#F39B7F", "Late" = "#E64B35")),
    show_legend = TRUE
  ),
  use_raster = TRUE
)
dev.off()

# ============================================================
# 4. 沿轨迹的细胞类型转换
# ============================================================
cat("\n--- 4. Cell Type Transitions ---\n")

# 把pseudotime分成5个bin，看每个bin中的细胞类型组成
plot_df_valid <- plot_df[valid_cells, ]
plot_df_valid <- plot_df_valid[ord, ]
plot_df_valid$pt_bin <- cut(plot_df_valid$pseudotime, breaks = 10, labels = FALSE)

comp <- plot_df_valid %>%
  group_by(pt_bin, cell_type) %>% summarise(n = n(), .groups = "drop") %>%
  group_by(pt_bin) %>% mutate(prop = n / sum(n))

p_comp <- ggplot(comp, aes(x = factor(pt_bin), y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors_fc) +
  theme_classic() +
  labs(title = "Cell Type Composition Along Pseudotime",
       x = "Pseudotime Bin (1=Early → 10=Late)", y = "Proportion", fill = "Cell Type") +
  theme(legend.text = element_text(size = 8))

pdf(file.path(res_dir, "figures/02_trajectory/celltype_composition_pseudotime.pdf"), width = 12, height = 6)
print(p_comp)
dev.off()

# ============================================================
# 5. 分化潜能分析（CytoTRACE-like）
# ============================================================
cat("\n--- 5. Differentiation Potential ---\n")

# 基于基因复杂度的分化潜能估算
gene_counts <- colSums(GetAssayData(obj, layer = "counts") > 0)
obj@meta.data$diff_potential <- gene_counts / max(gene_counts)

p_diff <- ggplot(obj@meta.data, aes(x = reorder(cell_type, -diff_potential), y = diff_potential, fill = cell_type)) +
  geom_boxplot(outlier.size = 0.2, show.legend = FALSE) +
  theme_classic() +
  labs(title = "Differentiation Potential by Cell Type",
       x = "", y = "Gene Complexity Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  scale_fill_manual(values = c(colors_fc,
    "Proliferating_G2M"="#66A61E", "Proliferating_M"="#A6D854", "Proliferating_S"="#B3DE69",
    "Male_Inflammatory"="#BC80BD", "Smooth_Muscle"="#8DD3C7",
    "Macrophage"="#4DBBD5", "Endothelial"="#00A087", "Stem_Progenitor"="#3C5488"))

pdf(file.path(res_dir, "figures/02_trajectory/differentiation_potential.pdf"), width = 12, height = 6)
print(p_diff)
dev.off()

# 保存
saveRDS(obj, file.path(data_dir, "processed/meniscus_annotated.rds"))

cat("\n=== Trajectory Analysis Complete ===\n")
cat("Figures saved to results/figures/02_trajectory/\n")
cat("  - slingshot_trajectory.pdf\n")
cat("  - pseudotime_umap.pdf\n")
cat("  - pseudotime_condition_density.pdf\n")
cat("  - gene_dynamics_heatmap.pdf\n")
cat("  - celltype_composition_pseudotime.pdf\n")
cat("  - differentiation_potential.pdf\n")
