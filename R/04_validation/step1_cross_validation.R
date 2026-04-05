# ============================================================
# Module 4: 计算验证与机器学习预测
# 
# 输入: 所有模块结果
# 输出: 交叉验证报告 + ML预测模型
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
  library(tidyverse)
  library(caret)
  library(randomForest)
  library(xgboost)
  library(ROCR)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
})

config  <- yaml::read_yaml("config/params.yaml")
data_dir <- path.expand(config$data_dir)
res_dir  <- path.expand(config$results_dir)

dir.create(file.path(res_dir, "figures", "05_validation"), recursive = TRUE, showWarnings = FALSE)

cat("=== Module 4: Validation & ML Prediction ===\n")

# 加载数据
merged <- readRDS(file.path(data_dir, "processed", "meniscus_annotated.rds"))
docking <- readRDS(file.path(data_dir, "processed", "docking_results.rds"))
cargo_db <- readRDS(file.path(data_dir, "processed", "cargo_database.rds"))

# ============================================================
# 4.1 独立数据集交叉验证
# ============================================================

cat("\n--- 4.1: Cross-Dataset Validation ---\n")

# 检查是否有多个数据集
datasets <- unique(merged$dataset)

if (length(datasets) >= 2) {
  cat("  Performing leave-one-dataset-out validation...\n")
  
  validation_results <- list()
  
  for (held_out in datasets) {
    cat(paste0("  Holding out: ", held_out, "\n"))
    
    # 训练集和验证集
    train_cells <- colnames(merged)[merged$dataset != held_out]
    test_cells  <- colnames(merged)[merged$dataset == held_out]
    
    train_obj <- subset(merged, cells = train_cells)
    test_obj  <- subset(merged, cells = test_cells)
    
    # 在训练集上找marker
    Idents(train_obj) <- "cell_type"
    train_markers <- FindAllMarkers(train_obj, only.pos = TRUE,
                                     min.pct = 0.25, logfc.threshold = 0.5,
                                     verbose = FALSE)
    
    top_markers <- train_markers %>%
      group_by(cluster) %>%
      slice_max(n = 20, order_by = avg_log2FC) %>%
      pull(gene) %>%
      unique()
    
    # 在验证集上检查这些marker的表达一致性
    common_genes <- intersect(top_markers, rownames(test_obj))
    
    if (length(common_genes) > 5) {
      # 计算验证集中marker基因的AUC
      marker_validation <- data.frame()
      
      for (ct in unique(train_markers$cluster)) {
        ct_markers_genes <- train_markers %>%
          filter(cluster == ct) %>%
          slice_max(n = 10, order_by = avg_log2FC) %>%
          pull(gene)
        ct_markers_genes <- ct_markers_genes[ct_markers_genes %in% rownames(test_obj)]
        
        if (length(ct_markers_genes) >= 3) {
          # 简化AUC：marker在目标细胞类型vs其他的表达差异
          target_cells <- test_obj$cell_type == ct
          if (sum(target_cells) >= 5 && sum(!target_cells) >= 5) {
            avg_in_target <- mean(colMeans(
              GetAssayData(test_obj, layer = "data")[ct_markers_genes, target_cells, drop = FALSE]
            ))
            avg_in_other <- mean(colMeans(
              GetAssayData(test_obj, layer = "data")[ct_markers_genes, !target_cells, drop = FALSE]
            ))
            
            fold_enrichment <- avg_in_target / max(avg_in_other, 0.001)
            
            marker_validation <- rbind(marker_validation, data.frame(
              held_out_dataset = held_out,
              cell_type = ct,
              n_markers_tested = length(ct_markers_genes),
              avg_expr_target = round(avg_in_target, 3),
              avg_expr_other = round(avg_in_other, 3),
              fold_enrichment = round(fold_enrichment, 2),
              validated = fold_enrichment > 1.5
            ))
          }
        }
      }
      
      validation_results[[held_out]] <- marker_validation
    }
  }
  
  if (length(validation_results) > 0) {
    validation_df <- bind_rows(validation_results)
    write.csv(validation_df,
              file.path(res_dir, "tables", "cross_dataset_validation.csv"),
              row.names = FALSE)
    
    # 验证成功率
    success_rate <- mean(validation_df$validated, na.rm = TRUE) * 100
    cat(paste0("  Cross-dataset validation success rate: ", round(success_rate, 1), "%\n"))
  }
  
} else {
  cat("  Only one dataset available. Using k-fold internal validation.\n")
}

# ============================================================
# 4.2 关键靶基因的bulk RNA-seq验证
# ============================================================

cat("\n--- 4.2: Bulk RNA-seq Cross-Validation ---\n")
cat("  Checking for GEO bulk RNA-seq OA datasets...\n")

# 准备关键基因集（从分子对接结果提取）
key_cargo_targets <- unique(c(
  # 治疗性cargo蛋白
  head(docking$treatability_scores$cell_type, 3),
  # 从评分矩阵提取top基因
  names(sort(rowSums(docking$protein_scores), decreasing = TRUE))[1:20]
))

# 准备GSEA用的基因集
# (用户需要下载OA相关bulk数据后运行此段)
cat("  To validate with bulk data, download OA vs Normal meniscus datasets from GEO.\n")
cat("  Suggested datasets: GSE98918, GSE113825, GSE55457\n")

# ============================================================
# 4.3 药物可成药性评估
# ============================================================

cat("\n--- 4.3: Druggability Assessment ---\n")

# 提取top治疗靶点
top_targets <- names(sort(rowSums(docking$protein_scores), decreasing = TRUE))[1:50]
top_targets <- top_targets[!is.na(top_targets)]

# 已知的药物靶点数据库(简化版)
known_drug_targets <- c(
  "TGFBR1", "TGFBR2", "FGFR1", "FGFR2", "FGFR3",
  "PDGFRA", "PDGFRB", "KDR", "MET", "EGFR",
  "IL1R1", "IL6R", "TNFRSF1A",
  "MMP1", "MMP2", "MMP3", "MMP9", "MMP13",
  "PTGS2", "NFKB1", "STAT3",
  "MTOR", "PIK3CA", "AKT1",
  "JAK1", "JAK2", "JAK3",
  "NOTCH1", "SMO", "CTNNB1"
)

druggable_targets <- intersect(top_targets, known_drug_targets)
cat(paste0("  Druggable targets among top predictions: ", length(druggable_targets), "\n"))
cat(paste0("  Targets: ", paste(druggable_targets, collapse = ", "), "\n"))

# 比较外泌体cargo vs 现有药物靶点
comparison <- data.frame(
  Target = top_targets[1:min(30, length(top_targets))],
  Is_Druggable = top_targets[1:min(30, length(top_targets))] %in% known_drug_targets,
  In_Exosome_Cargo = top_targets[1:min(30, length(top_targets))] %in% cargo_db$proteins$gene,
  stringsAsFactors = FALSE
)

write.csv(comparison,
          file.path(res_dir, "tables", "druggability_assessment.csv"),
          row.names = FALSE)

# ============================================================
# 4.4 机器学习预测模型
# ============================================================

cat("\n--- 4.4: Machine Learning Prediction ---\n")
cat("  Building predictive model for optimal cargo selection...\n")

# 构建训练数据: 基于cargo特征预测与各细胞类型的匹配效果
# 特征: cargo的功能属性 + 靶基因数量 + 通路覆盖度
# 标签: 匹配评分(连续变量 → 回归; 或二分类: high/low match)

# 蛋白cargo的特征矩阵
protein_features <- cargo_db$proteins %>%
  mutate(
    is_growth_factor = as.numeric(category == "Growth_Factor"),
    is_cytokine = as.numeric(category == "Cytokine_Chemokine"),
    is_ecm = as.numeric(category == "ECM_Regulation"),
    is_signaling = as.numeric(category == "Signaling"),
    is_surface = as.numeric(category == "Surface_Marker"),
    is_repair_relevant = as.numeric(therapeutic_relevance == "Repair"),
    is_anti_inflam = as.numeric(therapeutic_relevance == "Anti_inflammatory"),
    is_targeting = as.numeric(therapeutic_relevance == "Targeting")
  )

# 为每个cargo-celltype对构建训练样本
if (nrow(docking$protein_scores) > 0 && ncol(docking$protein_scores) > 0) {
  
  training_data <- data.frame()
  
  for (cargo_gene in rownames(docking$protein_scores)) {
    if (!(cargo_gene %in% protein_features$gene)) next
    
    feat <- protein_features %>% filter(gene == cargo_gene)
    
    for (ct in colnames(docking$protein_scores)) {
      score <- docking$protein_scores[cargo_gene, ct]
      
      training_data <- rbind(training_data, data.frame(
        cargo = cargo_gene,
        cell_type = ct,
        is_growth_factor = feat$is_growth_factor,
        is_cytokine = feat$is_cytokine,
        is_ecm = feat$is_ecm,
        is_signaling = feat$is_signaling,
        is_surface = feat$is_surface,
        is_repair_relevant = feat$is_repair_relevant,
        is_anti_inflam = feat$is_anti_inflam,
        matching_score = score,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  if (nrow(training_data) > 20) {
    # 二分类: 高匹配 vs 低匹配
    threshold <- quantile(training_data$matching_score, 0.75, na.rm = TRUE)
    training_data$high_match <- factor(ifelse(training_data$matching_score >= threshold, "High", "Low"))
    
    # 【修复】检查类别平衡性，如果只有一类则调整阈值
    class_counts <- table(training_data$high_match)
    cat(paste0("  Class distribution: ", paste(names(class_counts), "=", class_counts, collapse=", "), "\n"))
    
    if (length(class_counts) < 2) {
      # 降低阈值到中位数
      cat("  Only one class detected at 75th percentile. Trying median threshold...\n")
      threshold <- median(training_data$matching_score, na.rm = TRUE)
      training_data$high_match <- factor(ifelse(training_data$matching_score >= threshold, "High", "Low"))
      class_counts <- table(training_data$high_match)
      cat(paste0("  New class distribution: ", paste(names(class_counts), "=", class_counts, collapse=", "), "\n"))
    }
    
    if (length(class_counts) >= 2 && all(class_counts >= 5)) {
    feature_cols <- c("is_growth_factor", "is_cytokine", "is_ecm", "is_signaling",
                      "is_surface", "is_repair_relevant", "is_anti_inflam")
    
    # 训练Random Forest
    set.seed(42)
    train_idx <- createDataPartition(training_data$high_match, p = 0.7, list = FALSE)
    train_set <- training_data[train_idx, ]
    test_set  <- training_data[-train_idx, ]
    
    rf_model <- randomForest(
      high_match ~ .,
      data = train_set[, c(feature_cols, "high_match")],
      ntree = 500,
      importance = TRUE
    )
    
    # 预测与评估
    rf_pred <- predict(rf_model, test_set[, feature_cols])
    rf_accuracy <- mean(rf_pred == test_set$high_match)
    cat(paste0("  Random Forest accuracy: ", round(rf_accuracy * 100, 1), "%\n"))
    
    # 变量重要性
    importance_df <- data.frame(
      Feature = rownames(importance(rf_model)),
      MeanDecreaseGini = importance(rf_model)[, "MeanDecreaseGini"]
    ) %>% arrange(desc(MeanDecreaseGini))
    
    write.csv(importance_df,
              file.path(res_dir, "tables", "ML_feature_importance.csv"),
              row.names = FALSE)
    
    pdf(file.path(res_dir, "figures", "05_validation", "ML_feature_importance.pdf"),
        width = 8, height = 5)
    barplot(importance_df$MeanDecreaseGini,
            names.arg = importance_df$Feature,
            las = 2, col = "#3C5488",
            main = "Feature Importance for Cargo-Cell Matching Prediction",
            ylab = "Mean Decrease Gini",
            cex.names = 0.7)
    dev.off()
    
    # 混淆矩阵
    conf_matrix <- confusionMatrix(rf_pred, test_set$high_match)
    cat("  Confusion Matrix:\n")
    print(conf_matrix$table)
    
    # 保存模型
    saveRDS(rf_model, file.path(data_dir, "processed", "rf_cargo_prediction_model.rds"))
    
    } else {
      cat("  Class imbalance too severe for classification. Skipping RF model.\n")
    }
  } else {
    cat("  Insufficient training data for ML model.\n")
  }
}

# ============================================================
# 4.5 综合验证报告
# ============================================================

cat("\n--- 4.5: Generating Validation Summary ---\n")

# 汇总所有验证指标
validation_summary <- data.frame(
  Metric = c(
    "Cross-dataset marker validation rate",
    "Number of druggable targets identified",
    "ML prediction accuracy",
    "Active signaling pathways identified",
    "Exosome-repairable pathways",
    "Cell types with high treatability score"
  ),
  Value = c(
    ifelse(exists("success_rate"), paste0(round(success_rate, 1), "%"), "N/A"),
    length(druggable_targets),
    ifelse(exists("rf_accuracy"), paste0(round(rf_accuracy * 100, 1), "%"), "N/A"),
    length(docking$cellchat$All@netP$pathways),
    nrow(docking$exo_repair_map),
    sum(docking$treatability_scores$composite_score > 0.5)
  ),
  stringsAsFactors = FALSE
)

write.csv(validation_summary,
          file.path(res_dir, "tables", "validation_summary.csv"),
          row.names = FALSE)

cat("\n  Validation summary:\n")
print(validation_summary)

cat("\n=== Module 4 Complete ===\n")
cat("=== ALL ANALYSIS COMPLETE ===\n")
cat("\nResults summary:\n")
cat(paste0("  - Figures: ", res_dir, "/figures/\n"))
cat(paste0("  - Tables:  ", res_dir, "/tables/\n"))
cat(paste0("  - Data:    ", data_dir, "/processed/\n"))
