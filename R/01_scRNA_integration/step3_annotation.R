# ============================================================
# Step 3: 细胞类型注释与差异分析
# 
# 输入: 整合后Seurat对象
# 输出: 注释完成的对象 + DEG列表 + 可视化
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(yaml)
  library(tidyverse)
  library(patchwork)
  library(ComplexHeatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggsci)
})

config  <- yaml::read_yaml("config/params.yaml")
data_dir <- path.expand(config$data_dir)
res_dir  <- path.expand(config$results_dir)
markers  <- config$scRNA$cell_markers

dir.create(file.path(res_dir, "figures", "01_integration"), recursive = TRUE, showWarnings = FALSE)

cat("=== Step 3: Cell Type Annotation ===\n")

# ============================================================
# 3.1 加载数据
# ============================================================

merged <- readRDS(file.path(data_dir, "processed", "meniscus_integrated.rds"))
cat(paste0("  Loaded: ", ncol(merged), " cells, ", nlevels(Idents(merged)), " clusters\n"))

# ============================================================
# 3.2 标记基因鉴定
# ============================================================

cat("  Finding cluster markers...\n")

all_markers <- FindAllMarkers(merged,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.5,
  test.use = "wilcox",
  verbose = FALSE
)

# 保存完整marker表
write.csv(all_markers, file.path(res_dir, "tables", "cluster_markers_all.csv"),
          row.names = FALSE)

# Top markers per cluster
top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  ungroup()

write.csv(top_markers, file.path(res_dir, "tables", "cluster_markers_top10.csv"),
          row.names = FALSE)

# Top5 marker heatmap
top5 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene) %>%
  unique()

pdf(file.path(res_dir, "figures", "01_integration", "marker_heatmap.pdf"),
    width = 16, height = 12)
DoHeatmap(subset(merged, downsample = 100), features = top5, size = 3) +
  scale_fill_viridis_c()
dev.off()

# ============================================================
# 3.3 已知标记基因可视化
# ============================================================

cat("  Visualizing known markers...\n")

# 提取所有已知marker基因
known_markers <- unlist(markers, use.names = FALSE) %>% unique()
known_markers <- known_markers[known_markers %in% rownames(merged)]

# DotPlot
pdf(file.path(res_dir, "figures", "01_integration", "known_markers_dotplot.pdf"),
    width = 16, height = 8)
DotPlot(merged, features = known_markers, dot.scale = 6) +
  RotatedAxis() +
  scale_color_viridis_c() +
  theme(axis.text.x = element_text(size = 8))
dev.off()

# Feature plots for key markers
key_genes <- c("COL1A1", "COL2A1", "SOX9", "ACAN", "PRG4",
               "MCAM", "CD318", "MMP13", "PECAM1", "CD68",
               "PDGFRB", "ADAMTS5", "TGFB1", "IL1B")
key_genes <- key_genes[key_genes %in% rownames(merged)]

pdf(file.path(res_dir, "figures", "01_integration", "feature_plots_key_genes.pdf"),
    width = 20, height = 16)
FeaturePlot(merged, features = key_genes, ncol = 4, 
            cols = c("lightgrey", "red"), order = TRUE)
dev.off()

# ============================================================
# 3.4 自动注释 (SingleR)
# ============================================================

cat("  Running SingleR auto-annotation...\n")

# 使用 HumanPrimaryCellAtlas 作为参考
ref <- HumanPrimaryCellAtlasData()

# SingleR预测
sce <- as.SingleCellExperiment(merged)
singler_pred <- SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main,
  de.method = "wilcox"
)

merged$singler_labels <- singler_pred$labels

# ============================================================
# 3.5 综合注释（手动 + 自动）
# ============================================================

cat("  Performing manual annotation...\n")

# 基于marker基因的手动注释函数
annotate_cluster <- function(obj, cluster_id, marker_list) {
  cluster_markers_df <- all_markers %>% filter(cluster == cluster_id)
  top_genes <- cluster_markers_df$gene[1:50]
  
  scores <- sapply(names(marker_list), function(ct) {
    ct_markers <- marker_list[[ct]]
    overlap <- length(intersect(top_genes, ct_markers))
    avg_expr <- tryCatch({
      cells <- WhichCells(obj, idents = cluster_id)
      mean(colMeans(GetAssayData(obj, slot = "data")[ct_markers[ct_markers %in% rownames(obj)], cells, drop = FALSE]))
    }, error = function(e) 0)
    return(overlap * 0.4 + avg_expr * 0.6)
  })
  
  return(names(which.max(scores)))
}

# 对每个cluster执行注释
cluster_ids <- levels(Idents(merged))
annotations <- sapply(cluster_ids, function(cl) {
  tryCatch(
    annotate_cluster(merged, cl, markers),
    error = function(e) "Unknown"
  )
})

# 应用注释（添加编号以区分同类型的不同亚群）
annotation_counts <- table(annotations)
final_labels <- annotations
for (ct in names(annotation_counts[annotation_counts > 1])) {
  idx <- which(annotations == ct)
  for (i in seq_along(idx)) {
    final_labels[idx[i]] <- paste0(ct, "_", i)
  }
}

names(final_labels) <- cluster_ids
merged$cell_type <- final_labels[as.character(Idents(merged))]

# 如果自动注释有更好的结果，优先使用
# （此处保留手动注释结果，可根据实际情况调整）

cat("  Annotation results:\n")
print(table(merged$cell_type))

# ============================================================
# 3.6 注释后可视化
# ============================================================

cat("  Generating annotation plots...\n")

# 使用 NPG 配色
n_types <- length(unique(merged$cell_type))
colors <- pal_npg("nrc")(min(n_types, 10))
if (n_types > 10) {
  colors <- c(colors, pal_d3("category20")(n_types - 10))
}

# UMAP - cell types
p_ct <- DimPlot(merged, group.by = "cell_type", label = TRUE, repel = TRUE,
                label.size = 3.5, cols = colors) +
  ggtitle("Cell Type Annotation") +
  theme(legend.text = element_text(size = 8))

pdf(file.path(res_dir, "figures", "01_integration", "UMAP_cell_types.pdf"),
    width = 12, height = 8)
print(p_ct)
dev.off()

# 如果有条件信息(正常 vs 退变)，按条件分面
if ("condition" %in% colnames(merged@meta.data)) {
  p_split <- DimPlot(merged, group.by = "cell_type", split.by = "condition",
                     label = TRUE, repel = TRUE, cols = colors) +
    ggtitle("Cell Types by Condition")
  
  pdf(file.path(res_dir, "figures", "01_integration", "UMAP_split_condition.pdf"),
      width = 20, height = 8)
  print(p_split)
  dev.off()
}

# ============================================================
# 3.7 条件间差异表达分析
# ============================================================

cat("  Running differential expression between conditions...\n")

# 尝试识别条件信息
# (根据数据集的实际metadata列名调整)
condition_col <- NULL
for (col_name in c("condition", "group", "status", "disease", "type")) {
  if (col_name %in% colnames(merged@meta.data)) {
    condition_col <- col_name
    break
  }
}

if (!is.null(condition_col)) {
  conditions <- unique(merged@meta.data[[condition_col]])
  cat(paste0("  Found condition column: ", condition_col, " (", 
             paste(conditions, collapse = ", "), ")\n"))
  
  # 按细胞类型进行条件间差异分析
  deg_results <- list()
  for (ct in unique(merged$cell_type)) {
    ct_cells <- subset(merged, cell_type == ct)
    if (length(unique(ct_cells@meta.data[[condition_col]])) < 2) next
    if (min(table(ct_cells@meta.data[[condition_col]])) < 10) next
    
    tryCatch({
      Idents(ct_cells) <- condition_col
      degs <- FindMarkers(ct_cells,
        ident.1 = conditions[2],  # 退变/损伤
        ident.2 = conditions[1],  # 正常
        test.use = "wilcox",
        min.pct = 0.1,
        logfc.threshold = 0.25
      )
      degs$gene <- rownames(degs)
      degs$cell_type <- ct
      deg_results[[ct]] <- degs
    }, error = function(e) NULL)
  }
  
  if (length(deg_results) > 0) {
    deg_combined <- bind_rows(deg_results)
    write.csv(deg_combined,
              file.path(res_dir, "tables", "DEGs_by_condition_celltype.csv"),
              row.names = FALSE)
    cat(paste0("  Found DEGs in ", length(deg_results), " cell types\n"))
  }
  
} else {
  cat("  No condition column found. Skipping differential analysis.\n")
  cat("  You may need to add condition info manually based on sample metadata.\n")
}

# ============================================================
# 3.8 功能富集分析
# ============================================================

cat("  Running pathway enrichment...\n")

# 对每个cluster的marker基因做GO/KEGG富集
enrich_results <- list()
for (cl in unique(all_markers$cluster)) {
  cl_genes <- all_markers %>%
    filter(cluster == cl, p_val_adj < 0.05, avg_log2FC > 0.5) %>%
    pull(gene)
  
  if (length(cl_genes) < 10) next
  
  # 转换为Entrez ID
  gene_ids <- bitr(cl_genes, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
  
  if (nrow(gene_ids) < 5) next
  
  # GO enrichment
  ego <- enrichGO(
    gene = gene_ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  
  # KEGG enrichment
  ekegg <- enrichKEGG(
    gene = gene_ids$ENTREZID,
    organism = "hsa",
    pvalueCutoff = 0.05
  )
  
  enrich_results[[paste0("cluster_", cl)]] <- list(GO = ego, KEGG = ekegg)
}

# 保存富集结果
for (cl_name in names(enrich_results)) {
  if (!is.null(enrich_results[[cl_name]]$GO) && nrow(enrich_results[[cl_name]]$GO) > 0) {
    write.csv(as.data.frame(enrich_results[[cl_name]]$GO),
              file.path(res_dir, "tables", paste0("GO_enrichment_", cl_name, ".csv")),
              row.names = FALSE)
  }
}

# ============================================================
# 3.9 保存
# ============================================================

saveRDS(merged, file.path(data_dir, "processed", "meniscus_annotated.rds"))

cat("\n  Annotated object saved.\n")
cat("  Next step: Rscript scripts/01_scRNA_integration/step4_trajectory_GRN.R\n")
cat("=== Step 3 Complete ===\n")
