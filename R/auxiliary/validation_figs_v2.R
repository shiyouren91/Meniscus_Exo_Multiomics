# ============================================================
# 验证模块补充图表
# 上传到 AutoDL 后执行:
#   cd /root/autodl-tmp/meniscus_exosome_project
#   Rscript validation_figs_v2.R 2>&1 | tee logs/validation_figs.log
# ============================================================

options(conflicts.policy = list(warn = FALSE))
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(yaml)
  library(ComplexHeatmap); library(circlize)
})

config <- yaml::read_yaml("config/params.yaml")
data_dir <- config$data_dir; res_dir <- config$results_dir
dir.create(file.path(res_dir, "figures/05_validation"), recursive = TRUE, showWarnings = FALSE)

cat("=== Generating Validation Figures ===\n")

treat <- read.csv(file.path(res_dir, "tables/cell_type_treatability_scores.csv"))
treat$cell_type <- factor(treat$cell_type, levels = treat$cell_type[order(treat$composite_score)])

# 1. 可治疗性柱状图
cat("  1. Treatability barplot...\n")
pdf(file.path(res_dir, "figures/05_validation/treatability_barplot.pdf"), width = 10, height = 7)
print(ggplot(treat, aes(x = cell_type, y = composite_score)) +
  geom_bar(stat = "identity", aes(fill = composite_score), width = 0.7) +
  scale_fill_gradient2(low = "#4DBBD5", mid = "#F0E442", high = "#E64B35", midpoint = 0.45, guide = "none") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  coord_flip() + theme_classic() +
  labs(title = "Cell Type Treatability by MSC Exosomes", x = "", y = "Composite Treatability Score"))
dev.off()

# 2. 分解热图
cat("  2. Treatability decomposition...\n")
treat_mat <- as.matrix(treat[order(-treat$composite_score), c("receptor_richness","protein_cargo_match","mirna_target_expr","communication_involvement")])
rownames(treat_mat) <- treat$cell_type[order(-treat$composite_score)]
colnames(treat_mat) <- c("Receptor Richness","Protein Match","miRNA Target","Communication")

pdf(file.path(res_dir, "figures/05_validation/treatability_decomposition.pdf"), width = 10, height = 8)
Heatmap(treat_mat, name = "Score",
  col = colorRamp2(c(0, 0.5, 1), c("#F7F7F7", "#F39B7F", "#E64B35")),
  cluster_rows = FALSE, cluster_columns = FALSE,
  column_title = "Treatability Score Decomposition",
  cell_fun = function(j, i, x, y, w, h, fill) {
    v <- treat_mat[i, j]
    if (v > 0.01) grid.text(sprintf("%.2f", v), x, y, gp = gpar(fontsize = 8, col = if(v > 0.6) "white" else "black"))
  })
dev.off()

# 3. 差异受体火山图
cat("  3. Receptor volcano...\n")
de_rec <- read.csv(file.path(res_dir, "tables/differentially_expressed_receptors.csv"))
de_rec$sig <- case_when(
  de_rec$p_val_adj < 0.05 & de_rec$avg_log2FC > 0.5 ~ "Up",
  de_rec$p_val_adj < 0.05 & de_rec$avg_log2FC < -0.5 ~ "Down",
  TRUE ~ "NS")
de_rec$nlp <- pmin(-log10(de_rec$p_val + 1e-300), 50)

pdf(file.path(res_dir, "figures/05_validation/receptor_volcano.pdf"), width = 10, height = 8)
print(ggplot(de_rec, aes(avg_log2FC, nlp, color = sig)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = c(Up = "#E64B35", Down = "#4DBBD5", NS = "#CCCCCC")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  theme_classic() + labs(title = "Differential Receptors: Degenerated vs Normal",
                         x = "Log2FC", y = "-Log10(P)", color = ""))
dev.off()

cat("\n=== Validation Figures Complete ===\n")
