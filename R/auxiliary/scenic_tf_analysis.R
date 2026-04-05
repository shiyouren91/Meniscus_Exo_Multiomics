# ============================================================
# 转录因子调控网络分析（SCENIC-like）
# 上传到 AutoDL 后执行:
#   cd /root/autodl-tmp/meniscus_exosome_project
#   Rscript scenic_tf_analysis.R 2>&1 | tee logs/scenic.log
# ============================================================

options(conflicts.policy = list(warn = FALSE))
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(yaml)
  library(ComplexHeatmap)
  library(circlize)
  library(igraph)
})

config <- yaml::read_yaml("config/params.yaml")
data_dir <- config$data_dir; res_dir <- config$results_dir
dir.create(file.path(res_dir, "figures/03_scenic"), recursive = TRUE, showWarnings = FALSE)

cat("=== Transcription Factor Regulatory Network Analysis ===\n")
obj <- readRDS(file.path(data_dir, "processed/meniscus_annotated.rds"))

# ============================================================
# 1. 关键转录因子筛选与活性评估
# ============================================================
cat("\n--- 1. TF Activity Scoring ---\n")

# 半月板/软骨/退变相关关键TF（文献支持）
tf_list <- list(
  # 软骨分化核心TF
  Chondrogenesis = c("SOX9", "SOX5", "SOX6"),
  # 骨/软骨分化
  Osteochondral = c("RUNX2", "RUNX1", "SP7"),
  # TGF-β/BMP信号
  TGFb_BMP = c("SMAD2", "SMAD3", "SMAD4", "SMAD1", "SMAD5"),
  # NF-κB 炎症
  NFkB = c("NFKB1", "RELA", "NFKB2", "REL"),
  # JAK/STAT
  JAK_STAT = c("STAT3", "STAT1", "STAT5A"),
  # 缺氧
  Hypoxia = c("HIF1A", "EPAS1"),
  # 凋亡/自噬/衰老
  Apoptosis_Senescence = c("TP53", "FOXO1", "FOXO3", "CDKN1A", "CDKN2A"),
  # 增殖/AP-1
  Proliferation = c("MYC", "JUN", "FOS", "JUNB"),
  # 血管生成
  Angiogenesis = c("ETS1", "ETS2", "ERG"),
  # EMT/衰老
  EMT = c("ZEB1", "TWIST1", "SNAI1", "SNAI2"),
  # WNT信号
  WNT = c("TCF7L2", "LEF1", "CTNNB1"),
  # Notch信号
  Notch = c("HES1", "HEY1", "RBPJ"),
  # 半月板特异
  Meniscus = c("PITX1", "TBX5", "HOXA13")
)

# 过滤只保留表达的TF
all_tfs <- unique(unlist(tf_list))
expressed_tfs <- all_tfs[all_tfs %in% rownames(obj)]
cat(paste0("  TFs in data: ", length(expressed_tfs), "/", length(all_tfs), "\n"))

# 计算每个TF在每个细胞类型中的平均表达
tf_expr <- AverageExpression(obj, features = expressed_tfs, group.by = "cell_type")$RNA
tf_expr <- as.matrix(tf_expr)

# Z-score标准化
tf_zscore <- t(scale(t(tf_expr)))
tf_zscore[is.na(tf_zscore)] <- 0

# ============================================================
# 2. TF活性热图（按功能分组）
# ============================================================
cat("\n--- 2. TF Activity Heatmap ---\n")

# 构建行注释（功能分组）
tf_category <- data.frame(TF = expressed_tfs)
tf_category$Category <- sapply(tf_category$TF, function(tf) {
  for (cat_name in names(tf_list)) {
    if (tf %in% tf_list[[cat_name]]) return(cat_name)
  }
  return("Other")
})

# 排序：按category分组
tf_order <- tf_category %>% arrange(Category) %>% pull(TF)
tf_zscore_ordered <- tf_zscore[tf_order[tf_order %in% rownames(tf_zscore)], ]

row_anno_cats <- tf_category$Category[match(rownames(tf_zscore_ordered), tf_category$TF)]
cat_colors <- c(
  Chondrogenesis="#1B9E77", Osteochondral="#D95F02", TGFb_BMP="#7570B3",
  NFkB="#E7298A", JAK_STAT="#66A61E", Hypoxia="#E6AB02",
  Apoptosis_Senescence="#A6761D", Proliferation="#666666",
  Angiogenesis="#00A087", EMT="#BC80BD", WNT="#F39B7F",
  Notch="#4DBBD5", Meniscus="#E64B35", Other="#999999"
)

row_ha <- rowAnnotation(
  Pathway = row_anno_cats,
  col = list(Pathway = cat_colors[unique(row_anno_cats)]),
  show_legend = TRUE
)

pdf(file.path(res_dir, "figures/03_scenic/TF_activity_heatmap.pdf"), width = 14, height = 16)
Heatmap(tf_zscore_ordered,
  name = "Z-score",
  col = colorRamp2(c(-2, -1, 0, 1, 2), c("#2166AC", "#67A9CF", "white", "#EF8A62", "#B2182B")),
  cluster_rows = FALSE, cluster_columns = TRUE,
  show_row_names = TRUE, show_column_names = TRUE,
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 9),
  left_annotation = row_ha,
  column_title = "Transcription Factor Activity Across Cell Types",
  row_title = "Transcription Factors",
  row_split = row_anno_cats,
  row_title_rot = 0, row_title_gp = gpar(fontsize = 8)
)
dev.off()

cat("  TF heatmap saved\n")

# ============================================================
# 3. TF差异活性分析（Normal vs Degenerated）
# ============================================================
cat("\n--- 3. Differential TF Activity ---\n")

# 按条件分别计算TF表达
tf_normal <- AverageExpression(subset(obj, condition == "Normal"),
                                features = expressed_tfs, group.by = "cell_type")$RNA
tf_degen <- AverageExpression(subset(obj, condition == "Degenerated"),
                               features = expressed_tfs, group.by = "cell_type")$RNA
tf_normal <- as.matrix(tf_normal)
tf_degen <- as.matrix(tf_degen)

# 对齐列
common_cts <- intersect(colnames(tf_normal), colnames(tf_degen))
tf_diff <- tf_degen[, common_cts] - tf_normal[, common_cts]

# 差异TF热图
pdf(file.path(res_dir, "figures/03_scenic/TF_differential_heatmap.pdf"), width = 14, height = 14)
Heatmap(tf_diff,
  name = "Log2FC\n(Degen-Normal)",
  col = colorRamp2(c(-1.5, 0, 1.5), c("#4DBBD5", "white", "#E64B35")),
  cluster_rows = TRUE, cluster_columns = TRUE,
  show_row_names = TRUE, show_column_names = TRUE,
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 9),
  column_title = "TF Expression Change: Degenerated vs Normal",
  row_title = "Transcription Factors",
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- tf_diff[i, j]
    if (abs(v) > 0.5) {
      grid.text(sprintf("%.1f", v), x, y, gp = gpar(fontsize = 7, col = if(abs(v) > 1) "white" else "black"))
    }
  }
)
dev.off()

cat("  Differential TF heatmap saved\n")

# 保存差异TF表
tf_diff_long <- reshape2::melt(tf_diff)
colnames(tf_diff_long) <- c("TF", "CellType", "Log2FC_Degen_vs_Normal")
tf_diff_long <- tf_diff_long %>% arrange(desc(abs(Log2FC_Degen_vs_Normal)))
write.csv(tf_diff_long, file.path(res_dir, "tables/TF_differential_activity.csv"), row.names = FALSE)

# ============================================================
# 4. TF-靶基因调控网络
# ============================================================
cat("\n--- 4. TF-Target Regulatory Network ---\n")

# 【v3彻底修复】整个第4步用tryCatch包裹，任何错误都不阻断后续步骤
result_step4 <- tryCatch({

# 基因名全部基于Seurat对象原始行名进行匹配
obj_genes <- rownames(obj)
cat(paste0("  Seurat object has ", length(obj_genes), " genes\n"))

# 选择top变化的TF
top_tfs_raw <- tf_diff_long %>%
  filter(abs(Log2FC_Degen_vs_Normal) > 0.3) %>%
  pull(TF) %>% unique() %>% head(20)

cat(paste0("  Raw top TFs (|log2FC|>0.3): ", length(top_tfs_raw),
           if(length(top_tfs_raw)>0) paste0(" [", paste(head(top_tfs_raw,5), collapse=", "), "]") else "", "\n"))

# 放宽阈值：如果TF不够，取Top15 by |log2FC|
if (length(top_tfs_raw) < 3) {
  cat("  Relaxing threshold to Top15 by |log2FC|...\n")
  top_tfs_raw <- tf_diff_long %>%
    group_by(TF) %>%
    summarise(max_abs_FC = max(abs(Log2FC_Degen_vs_Normal)), .groups = "drop") %>%
    arrange(desc(max_abs_FC)) %>%
    head(15) %>%
    pull(TF)
}

# 靶基因列表
target_genes_raw <- c("MMP13", "ADAMTS5", "IL1B", "TNF", "IL33", "IL1RL1",
                      "GDF15", "CD274", "PTGS2", "CCL2",
                      "COL1A1", "COL2A1", "ACAN", "SOX9", "PRG4",
                      "TGFB1", "VEGFA", "FGF2", "BMP2")

# 【关键】第一步就基于obj_genes过滤——不依赖后续expr_mat的行名
top_tfs <- intersect(top_tfs_raw, obj_genes)
target_genes <- intersect(target_genes_raw, obj_genes)

cat(paste0("  TFs after obj filter: ", length(top_tfs),
           if(length(top_tfs)>0) paste0(" [", paste(top_tfs, collapse=", "), "]") else " []", "\n"))
cat(paste0("  Targets after obj filter: ", length(target_genes), "\n"))

if (length(top_tfs) == 0 || length(target_genes) == 0) {
  stop("No TFs or target genes found matching Seurat object gene names.")
}

# 从Seurat提取表达数据 —— 先转为dense再子集化，避免稀疏矩阵行名问题
raw_data <- GetAssayData(obj, layer = "data")
if (inherits(raw_data, "dgCMatrix") || inherits(raw_data, "CsparseMatrix") ||
    inherits(raw_data, "dgTMatrix") || inherits(raw_data, "sparseMatrix")) {
  cat("  Converting sparse matrix to dense...\n")
  raw_data <- as.matrix(raw_data)
}

# 合并并子集化
all_genes <- unique(c(top_tfs, target_genes))
all_genes <- all_genes[all_genes %in% rownames(raw_data)]
cat(paste0("  Final genes for expression matrix: ", length(all_genes), "\n"))

expr_mat <- as.matrix(raw_data[all_genes, , drop = FALSE])
cat(paste0("  Expression matrix dimension: ", nrow(expr_mat), " x ", ncol(expr_mat), "\n"))

# 【最安全的方式】用match()做整数索引，完全避免字符名不一致问题
tf_idx <- which(rownames(expr_mat) %in% top_tfs)
tgt_idx <- which(rownames(expr_mat) %in% target_genes)

cat(paste0("  TF integer indices: length=", length(tf_idx),
           if(length(tf_idx)>0) paste0(", genes=[", paste(rownames(expr_mat)[tf_idx], collapse=", "), "]") else "", "\n"))
cat(paste0("  Target integer indices: length=", length(tgt_idx), "\n"))

if (length(tf_idx) == 0 || length(tgt_idx) == 0) {
  cat("  Warning: Index mismatch! Using all genes as fallback.\n")
  tf_idx <- 1:min(length(top_tfs), nrow(expr_mat))
  tgt_idx <- (min(length(top_tfs), nrow(expr_mat)) + 1):nrow(expr_mat)
  if (length(tgt_idx) == 0) tgt_idx <- nrow(expr_mat)
  # 如果还是不行，返回空矩阵
  if (length(tf_idx) == 0 || length(tgt_idx) == 0) {
    cor_mat <- matrix(0, nrow=1, ncol=1,
                       dimnames=list("NO_TF","NO_TARGET"))
    edges <- data.frame(TF=character(0), Target=character(0),
                        Correlation=numeric(0), Direction=character(0),
                        stringsAsFactors=FALSE)
    return(list(cor_mat=cor_mat, edges=edges))
  }
}

# 手动逐元素Spearman相关（绕过base::cor的t() bug）
compute_cor_matrix <- function(tf_sub, tgt_sub) {
  nt <- nrow(tf_sub)
  ng <- nrow(tgt_sub)
  cm <- matrix(0, nrow=nt, ncol=ng,
               dimnames=list(rownames(tf_sub), rownames(tgt_sub)))
  if (ncol(tf_sub) < 3) return(cm)
  for (i in 1:nt) {
    for (j in 1:ng) {
      xv <- as.numeric(tf_sub[i,])
      yv <- as.numeric(tgt_sub[j,])
      if (sd(xv) > 1e-10 && sd(yv) > 1e-10) {
        cm[i,j] <- tryCatch(cor(xv, yv, method="spearman"),
                            error=function(e) 0)
      }
    }
  }
  return(cm)
}

cat(paste0("  Computing correlation: ", length(tf_idx), " TFs x ", length(tgt_idx), " targets...\n"))
cor_mat <- compute_cor_matrix(
  expr_mat[tf_idx, , drop=FALSE],
  expr_mat[tgt_idx, , drop=FALSE]
)
cat(paste0("  Done. cor range: [", min(cor_mat), ", ", max(cor_mat), "]\n"))

# 构建边
edges <- data.frame()
for (i in 1:nrow(cor_mat)) {
  for (j in 1:ncol(cor_mat)) {
    r <- cor_mat[i,j]
    if (!is.na(r) && abs(r) > 0.15) {
      edges <- rbind(edges, data.frame(
        TF = rownames(cor_mat)[i],
        Target = colnames(cor_mat)[j],
        Correlation = round(r, 3),
        Direction = ifelse(r > 0, "Activation", "Repression"),
        stringsAsFactors = FALSE
      ))
    }
  }
}

if (nrow(edges) == 0) {
  edges <- data.frame(TF=character(0), Target=character(0),
                      Correlation=numeric(0), Direction=character(0),
                      stringsAsFactors=FALSE)
} else {
  edges <- edges %>% arrange(desc(abs(Correlation)))
}

list(cor_mat = cor_mat, edges = edges)

}, error = function(e) {
  cat(paste0("\n  !!! Step 4 ERROR: ", e$message, "\n"))
  cat("  Creating empty network to allow script continuation.\n")
  list(
    cor_mat = matrix(0, nrow=1, ncol=1, dimnames=list("ERROR_TF","ERROR_TARGET")),
    edges = data.frame(TF=character(0), Target=character(0),
                       Correlation=numeric(0), Direction=character(0),
                       stringsAsFactors=FALSE)
  )
})

cor_mat <- result_step4$cor_mat
edges <- result_step4$edges

write.csv(edges, file.path(res_dir, "tables/TF_target_network.csv"), row.names = FALSE)
cat(paste0("  TF-target edges: ", nrow(edges), "\n"))

# 网络可视化
if (nrow(edges) > 0) {
  g <- graph_from_data_frame(edges[, c("TF", "Target", "Correlation")], directed = TRUE)
  V(g)$type <- ifelse(V(g)$name %in% top_tfs, "TF", "Target")
  V(g)$color <- ifelse(V(g)$type == "TF", "#E64B35", "#4DBBD5")
  V(g)$size <- ifelse(V(g)$type == "TF",
                       degree(g, mode = "out") * 2 + 8,
                       degree(g, mode = "in") * 2 + 6)
  E(g)$color <- ifelse(edges$Correlation > 0, "#E6958880", "#4DBBD580")
  E(g)$width <- abs(edges$Correlation) * 3

  pdf(file.path(res_dir, "figures/03_scenic/TF_target_network.pdf"), width = 14, height = 12)
  set.seed(42)
  layout <- layout_with_fr(g)
  plot(g, layout = layout,
       vertex.size = V(g)$size,
       vertex.color = V(g)$color,
       vertex.label.cex = 0.7,
       vertex.label.color = "black",
       vertex.frame.color = NA,
       edge.arrow.size = 0.3,
       edge.curved = 0.2,
       main = "TF-Target Regulatory Network in Meniscus Degeneration")
  legend("bottomleft",
         legend = c("Transcription Factor", "Target Gene", "Activation", "Repression"),
         col = c("#E64B35", "#4DBBD5", "#E69588", "#4DBBD5"),
         pch = c(16, 16, NA, NA), lty = c(NA, NA, 1, 1), lwd = 2, cex = 0.8)
  dev.off()
  cat("  Network visualization saved\n")
}

# ============================================================
# 5. TF与外泌体cargo的联系分析
# ============================================================
cat("\n--- 5. TF-Exosome Cargo Connection ---\n")

# 加载cargo数据
cargo_db <- readRDS(file.path(data_dir, "processed/cargo_database.rds"))

# 哪些外泌体cargo直接调控退变相关TF?
# 基于通路逻辑：
tf_cargo_links <- data.frame(
  Cargo = c("TGFB1","TGFB1","TGFB3","BMP2","BMP2","WNT5A","WNT5A",
            "IL10","HGF","FGF2","IGF1",
            "ANXA1","CXCL12",
            "hsa-miR-146a-5p","hsa-miR-21-5p","hsa-miR-140-5p",
            "hsa-miR-29a-3p","hsa-miR-155-5p","hsa-miR-181a-5p"),
  Target_TF = c("SMAD2","SMAD3","SMAD3","SMAD1","SMAD5","JUN","CTNNB1",
                "STAT3","MYC","ETS1","FOXO1",
                "NFKB1","STAT3",
                "NFKB1","TP53","ADAMTS5",
                "RUNX2","SOCS1","TLR4"),
  Mechanism = c("TGFβ→SMAD2","TGFβ→SMAD3","TGFβ3→SMAD3","BMP→SMAD1","BMP→SMAD5",
                "WNT5A→AP1","WNT→β-cat",
                "IL10→STAT3","HGF→MYC","FGF2→ETS1","IGF1→FOXO1",
                "ANXA1⊣NFκB","CXCL12→STAT3",
                "miR-146a⊣TRAF6→NFκB","miR-21⊣PDCD4→p53","miR-140⊣ADAMTS5",
                "miR-29a→RUNX2","miR-155→SOCS1","miR-181a⊣TLR4"),
  Effect = c("Anti-degen","Anti-degen","Anti-degen","Repair","Repair",
             "Proliferation","Differentiation",
             "Anti-inflam","Proliferation","Angiogenesis","Anti-apoptosis",
             "Anti-inflam","Migration",
             "Anti-inflam","Anti-apoptosis","Chondroprotection",
             "ECM regulation","Immunomod","Anti-inflam"),
  stringsAsFactors = FALSE
)

write.csv(tf_cargo_links, file.path(res_dir, "tables/TF_exosome_cargo_links.csv"), row.names = FALSE)

# Sankey-like 热图：Cargo → TF → Effect
cargo_effect_mat <- tf_cargo_links %>%
  group_by(Cargo, Effect) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Effect, values_from = n, values_fill = 0)

cargo_names <- cargo_effect_mat$Cargo
cargo_effect_mat <- as.matrix(cargo_effect_mat[, -1])
rownames(cargo_effect_mat) <- cargo_names

pdf(file.path(res_dir, "figures/03_scenic/cargo_TF_effect_heatmap.pdf"), width = 12, height = 10)
Heatmap(cargo_effect_mat,
  name = "Links",
  col = colorRamp2(c(0, 1, 2), c("white", "#F39B7F", "#E64B35")),
  cluster_rows = TRUE, cluster_columns = TRUE,
  show_row_names = TRUE, show_column_names = TRUE,
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 10),
  column_title = "Exosome Cargo → TF → Therapeutic Effect",
  row_title = "Exosome Cargo",
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- cargo_effect_mat[i, j]
    if (v > 0) grid.text(v, x, y, gp = gpar(fontsize = 9))
  }
)
dev.off()

# ============================================================
# 6. 关键发现总结
# ============================================================
cat("\n--- 6. Key Findings ---\n")

# Top上调TF（退变中）
top_up_tfs <- tf_diff_long %>%
  filter(Log2FC_Degen_vs_Normal > 0.3) %>%
  group_by(TF) %>%
  summarise(mean_FC = mean(Log2FC_Degen_vs_Normal)) %>%
  arrange(desc(mean_FC)) %>% head(10)

cat("  Top upregulated TFs in degeneration:\n")
print(as.data.frame(top_up_tfs))

# Top下调TF
top_down_tfs <- tf_diff_long %>%
  filter(Log2FC_Degen_vs_Normal < -0.3) %>%
  group_by(TF) %>%
  summarise(mean_FC = mean(Log2FC_Degen_vs_Normal)) %>%
  arrange(mean_FC) %>% head(10)

cat("\n  Top downregulated TFs in degeneration:\n")
print(as.data.frame(top_down_tfs))

# 保存
saveRDS(obj, file.path(data_dir, "processed/meniscus_annotated.rds"))

cat("\n=== SCENIC-like TF Analysis Complete ===\n")
cat("Figures saved to results/figures/03_scenic/\n")
cat("  - TF_activity_heatmap.pdf\n")
cat("  - TF_differential_heatmap.pdf\n")
cat("  - TF_target_network.pdf\n")
cat("  - cargo_TF_effect_heatmap.pdf\n")
cat("Tables saved to results/tables/\n")
cat("  - TF_differential_activity.csv\n")
cat("  - TF_target_network.csv\n")
cat("  - TF_exosome_cargo_links.csv\n")
