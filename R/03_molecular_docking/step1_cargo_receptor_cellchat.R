# ============================================================
# Module 3: 跨尺度分子对接（核心创新模块）
# Step 1+2+3 合并: 受体指纹 → Cargo-受体匹配 → CellChat外泌体通讯
#
# 输入: 注释后的Seurat对象 + Cargo数据库 + miRNA靶基因网络
# 输出: 外泌体-半月板跨尺度分子对接图谱
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(yaml)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(igraph)
  library(ggalluvial)
  library(patchwork)
  library(ggsci)
  library(openxlsx)
})

config  <- yaml::read_yaml("config/params.yaml")
data_dir <- path.expand(config$data_dir)
res_dir  <- path.expand(config$results_dir)
scoring_weights <- config$docking$scoring_weights

dir.create(file.path(res_dir, "figures", "04_docking"), recursive = TRUE, showWarnings = FALSE)

cat("=== Module 3: Cross-Scale Molecular Docking ===\n")

# 加载数据
merged   <- readRDS(file.path(data_dir, "processed", "meniscus_annotated.rds"))
cargo_db <- readRDS(file.path(data_dir, "processed", "cargo_database.rds"))
mirna_analysis <- readRDS(file.path(data_dir, "processed", "mirna_target_analysis.rds"))

# ============================================================
# PART 1: 受体端分析 — 细胞亚群特异性受体指纹
# ============================================================

cat("\n--- Part 1: Receptor Fingerprinting ---\n")

# 1.1 定义受体/靶基因集合
# 从CellChatDB提取所有受体
CellChatDB <- CellChatDB.human
all_receptors <- unique(c(
  unlist(strsplit(CellChatDB$interaction$receptor, "_")),
  # 额外的半月板相关受体
  "TGFBR1", "TGFBR2", "BMPR1A", "BMPR2", "ACVR1",
  "FGFR1", "FGFR2", "FGFR3", "IGFR1",
  "IL1R1", "IL6R", "TNFRSF1A", "TNFRSF1B",
  "ITGA5", "ITGB1", "ITGAV", "ITGB3",
  "CD44", "SDC1", "SDC4",
  "PDGFRA", "PDGFRB", "KDR", "FLT1",
  "MET", "EGFR", "ERBB2",
  "NOTCH1", "NOTCH2", "NOTCH3",
  "FZD1", "FZD2", "FZD5", "LRP5", "LRP6",
  "PTCH1", "SMO"
))

# 过滤只保留表达的受体
expressed_receptors <- all_receptors[all_receptors %in% rownames(merged)]
cat(paste0("  Expressed receptors: ", length(expressed_receptors), "\n"))

# 1.2 计算各细胞类型的受体表达谱
receptor_expr <- AverageExpression(merged, 
  features = expressed_receptors,
  group.by = "cell_type"
)$RNA

# 标准化 (Z-score per receptor)
receptor_zscore <- t(scale(t(as.matrix(receptor_expr))))
receptor_zscore[is.na(receptor_zscore)] <- 0

# 1.3 受体指纹热图
pdf(file.path(res_dir, "figures", "04_docking", "receptor_fingerprint_heatmap.pdf"),
    width = 14, height = 16)
Heatmap(
  receptor_zscore,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("#4DBBD5", "white", "#E64B35")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 7),
  column_title = "Cell Type-Specific Receptor Fingerprint",
  row_title = "Receptors"
)
dev.off()

# 1.4 识别损伤后差异表达的受体
# (如果有条件对比数据)
condition_col <- NULL
for (col_name in c("condition", "group", "status")) {
  if (col_name %in% colnames(merged@meta.data)) {
    condition_col <- col_name
    break
  }
}

deg_receptors <- NULL
if (!is.null(condition_col)) {
  cat("  Identifying differentially expressed receptors...\n")
  
  deg_receptors <- list()
  for (ct in unique(merged$cell_type)) {
    ct_cells <- subset(merged, cell_type == ct)
    conditions <- unique(ct_cells@meta.data[[condition_col]])
    if (length(conditions) < 2 || min(table(ct_cells@meta.data[[condition_col]])) < 10) next
    
    Idents(ct_cells) <- condition_col
    tryCatch({
      degs <- FindMarkers(ct_cells, ident.1 = conditions[2], ident.2 = conditions[1],
                          features = expressed_receptors, min.pct = 0.05,
                          logfc.threshold = 0.1)
      degs$receptor <- rownames(degs)
      degs$cell_type <- ct
      deg_receptors[[ct]] <- degs
    }, error = function(e) NULL)
  }
  
  if (length(deg_receptors) > 0) {
    deg_receptors_df <- bind_rows(deg_receptors)
    write.csv(deg_receptors_df,
              file.path(res_dir, "tables", "differentially_expressed_receptors.csv"),
              row.names = FALSE)
  }
}

# ============================================================
# PART 2: Cargo-受体匹配矩阵
# ============================================================

cat("\n--- Part 2: Cargo-Receptor Matching ---\n")

# 2.1 蛋白Cargo → 受体配对 (基于CellChatDB)
cat("  Matching protein cargo to receptors...\n")

cargo_proteins <- cargo_db$proteins$gene
cat(paste0("  Total cargo proteins: ", length(cargo_proteins), "\n"))
cat(paste0("  Cargo proteins: [", paste(head(cargo_proteins, 10), collapse=", "), 
           if(length(cargo_proteins)>10) "..." else "", "]\n"))

interaction_db <- CellChatDB$interaction
cat(paste0("  CellChatDB interactions total: ", nrow(interaction_db), "\n"))

# 【调试】展示CellChatDB ligand/receptor名称示例
if (nrow(interaction_db) > 0) {
  cat(paste0("  Sample ligands: [", paste(head(interaction_db$ligand, 5), collapse=", "), "]\n"))
  cat(paste0("  Sample receptors: [", paste(head(interaction_db$receptor, 5), collapse=", "), "]\n"))
}

# 找到cargo蛋白作为配体的interaction
# 【修复】使用更宽松的匹配策略：部分匹配 + 大小写不敏感
cargo_ligand_matches <- interaction_db %>%
  filter(sapply(strsplit(ligand, "_"), function(x) {
    any(tolower(x) %in% tolower(cargo_proteins))
  }))

cat(paste0("  Cargo-receptor pairs from CellChatDB: ", nrow(cargo_ligand_matches), "\n"))

if (nrow(cargo_ligand_matches) == 0) {
  cat("  !!! No matches found via exact name. Trying fuzzy matching...\n")
  
  # 【回退方案】手动构建已知的cargo-receptor相互作用网络（基于文献）
  # 这些是MSC外泌体cargo与其受体的已知配对，不依赖CellChatDB名称格式
  known_cargo_receptor_pairs <- data.frame(
    ligand = c(
      "TGFB1", "TGFB3", "TGFB2",
      "BMP2", "BMP7",
      "FGF2", "FGF1", "HGF",
      "IGF1", "VEGFA",
      "IL10", "CCL2", "CXCL12",
      "PDGFB",
      "WNT3A", "WNT5A"
    ),
    receptor = c(
      "TGFBR1_TGFBR2", "TGFBR1_TGFBR2", "TGFBR1_TGFBR2",
      "BMPR1A_BMPR2", "BMPR1A_BMPR2",
      "FGFR1_FGFR2", "FGFR1_FGFR2", "MET",
      "IGFR1", "KDR_FLT1",
      "IL10RA_IL10RB", "CCR2", "CXCR4",
      "PDGFRA_PDGFRB",
      "FZD1_FZD2_LRP5_LRP6", "FZD1_FZD2_LRP5_LRP6"
    ),
    pathway_name = c(
      "TGF-beta", "TGF-beta", "TGF-beta",
      "BMP", "BMP",
      "FGF", "FGF", "HGF",
      "IGF", "VEGF",
      "Inflammatory", "Inflammatory", "Inflammatory",
      "PDGF",
      "WNT", "WNT"
    ),
    stringsAsFactors = FALSE
  )
  
  # 只保留我们cargo库中有的蛋白
  known_cargo_receptor_pairs <- known_cargo_receptor_pairs %>%
    filter(ligand %in% cargo_proteins)
  
  cat(paste0("  Fallback: Using ", nrow(known_cargo_receptor_pairs), 
             " literature-curated cargo-receptor pairs\n"))
  
  if (nrow(known_cargo_receptor_pairs) > 0) {
    cargo_ligand_matches <- known_cargo_receptor_pairs
  }
} else {
  # 展示匹配到的具体pairs
  if (nrow(cargo_ligand_matches) <= 30) {
    for (i in 1:min(10, nrow(cargo_ligand_matches))) {
      cat(paste0("    Match ", i, ": ", cargo_ligand_matches$ligand[i], 
                 " -> ", cargo_ligand_matches$receptor[i],
                 " [", cargo_ligand_matches$pathway_name[i], "]\n"))
    }
  } else {
    cat(paste0("    Showing first 10 of ", nrow(cargo_ligand_matches), " matches\n"))
    for (i in 1:10) {
      cat(paste0("    Match ", i, ": ", cargo_ligand_matches$ligand[i],
                 " -> ", cargo_ligand_matches$receptor[i], "\n"))
    }
  }
}

# 2.2 miRNA Cargo → 靶基因 → 细胞类型分配
cat("  Mapping miRNA targets to cell types...\n")

mirna_targets <- mirna_analysis$high_conf_targets

# 检查靶基因在各细胞类型中的表达
target_genes <- unique(mirna_targets$Target)
target_genes <- target_genes[target_genes %in% rownames(merged)]

target_expr_by_ct <- AverageExpression(merged, 
  features = target_genes,
  group.by = "cell_type"
)$RNA

# 2.3 构建 Cargo × CellType 匹配评分矩阵
cat("  Building Cargo × CellType scoring matrix...\n")

cell_types <- unique(merged$cell_type)

# --- 蛋白Cargo评分 ---
protein_scores <- matrix(0, nrow = length(cargo_proteins), ncol = length(cell_types),
                          dimnames = list(cargo_proteins, cell_types))

# 【调试】检查receptor_expr的行名格式
cat(paste0("  receptor_expr dim: ", nrow(receptor_expr), " x ", ncol(receptor_expr), "\n"))
cat(paste0("  receptor_expr rownames sample: [", paste(head(rownames(receptor_expr), 5), collapse=", "), "]\n"))
cat(paste0("  expressed_receptors (", length(expressed_receptors), "): [",
           paste(head(expressed_receptors, 8), collapse=", "), 
           if(length(expressed_receptors)>8) "..." else "", "]\n"))

# 【关键修复】大小写不敏感的receptor匹配
# receptor_expr可能使用不同的大小写格式
receptor_expr_rownames_lower <- tolower(rownames(receptor_expr))
expressed_receptors_lower <- tolower(expressed_receptors)

n_matched_pairs <- 0
n_zero_signal_pairs <- 0

for (i in seq_len(nrow(cargo_ligand_matches))) {
  ligand_genes <- unlist(strsplit(as.character(cargo_ligand_matches$ligand[i]), "_"))
  receptor_genes <- unlist(strsplit(as.character(cargo_ligand_matches$receptor[i]), "_"))
  
  matched_cargo <- intersect(ligand_genes, cargo_proteins)
  
  # 【修复】先用精确匹配，再用大小写不敏感匹配
  matched_receptors <- intersect(receptor_genes, expressed_receptors)
  if (length(matched_receptors) == 0 && length(receptor_genes) > 0) {
    # 大小写不敏感回退
    lower_matches <- which(tolower(receptor_genes) %in% expressed_receptors_lower)
    if (length(lower_matches) > 0) {
      matched_receptors <- receptor_genes[lower_matches]
    }
  }
  
  if (length(matched_cargo) > 0 && length(matched_receptors) > 0) {
    # 【v5修复】对每个receptor做完整的安全解析链：
    #   精确匹配 → 大小写不敏感 → 复合名拆分 → 最终安全取值
    valid_receptor_names <- character(0)
    
    for (rec in matched_receptors) {
      rec_name <- NULL
      # Step A: 精确匹配
      if (rec %in% rownames(receptor_expr)) {
        rec_name <- rec
      } else {
        # Step B: 大小写不敏感匹配
        lower_rec <- tolower(rec)
        expr_lower <- tolower(rownames(receptor_expr))
        fuzzy_hits <- which(expr_lower == lower_rec)
        if (length(fuzzy_hits) > 0) {
          rec_name <- rownames(receptor_expr)[fuzzy_hits[1]]
          cat(paste0("    FIX case: ", rec, " -> ", rec_name, "\n"))
        } else {
          # Step C: 尝试从复合名中找子串匹配（如 TGFbR1_R2 -> TGFBR1）
          sub_parts <- unlist(strsplit(rec, "_"))
          for (sp in sub_parts) {
            if (sp %in% rownames(receptor_expr)) {
              rec_name <- sp
              break
            }
            sp_lower <- tolower(sp)
            sub_fuzzy <- which(expr_lower == sp_lower)
            if (length(sub_fuzzy) > 0) {
              rec_name <- rownames(receptor_expr)[sub_fuzzy[1]]
              cat(paste0("    FIX substring: ", sp, " -> ", rec_name, " (from ", rec, ")\n"))
              break
            }
          }
        }
      }
      if (!is.null(rec_name)) {
        valid_receptor_names <- c(valid_receptor_names, rec_name)
      } else {
        if (i <= 5) cat(paste0("    SKIP: receptor '", rec, "' not found in expression matrix\n"))
      }
    }
    
    if (length(valid_receptor_names) > 0) {
      for (ct in cell_types) {
        # 【v5修复】逐个取值再mean，避免向量下标越界
        signal_values <- numeric(0)
        for (rn in valid_receptor_names) {
          val <- tryCatch({
            v <- receptor_expr[rn, ct]
            if (length(v) > 0 && !is.na(v)) as.numeric(v) else NA
          }, error = function(e) { 
            NA 
          })
          if (!is.na(val)) signal_values <- c(signal_values, val)
        }
        
        receptor_signal <- if (length(signal_values) > 0) mean(signal_values) else NA
        
        if (!is.na(receptor_signal) && receptor_signal > 0) {
          for (cargo in matched_cargo) {
            protein_scores[cargo, ct] <- protein_scores[cargo, ct] + receptor_signal
          }
          n_matched_pairs <- n_matched_pairs + 1
        } else {
          n_zero_signal_pairs <- n_zero_signal_pairs + 1
        }
      }
    } else {
      cat(paste0("    WARN: Receptors found but not in receptor_expr: ",
                 paste(matched_receptors, collapse=", "), "\n"))
    }
  } else {
    if (length(matched_cargo) > 0 && i <= 5) {
      cat(paste0("    No receptors matched for cargo ligand: ", 
                 cargo_ligand_matches$ligand[i], " (receptors: ",
                 paste(receptor_genes, collapse=", "), ")\n"))
    }
  }
}

cat(paste0("  Scored pairs with signal: ", n_matched_pairs, ", zero-signal: ", n_zero_signal_pairs, "\n"))
cat(paste0("  Non-zero entries in protein_scores: ", sum(protein_scores > 0), "/", length(protein_scores), "\n"))
cat(paste0("  protein_scores range: [", min(protein_scores), ", ", max(protein_scores), "]\n"))

# --- miRNA Cargo评分 ---
mirna_ids <- unique(mirna_targets$miRNA)
mirna_scores <- matrix(0, nrow = length(mirna_ids), ncol = length(cell_types),
                        dimnames = list(mirna_ids, cell_types))

for (mir in mirna_ids) {
  targets <- mirna_targets %>% filter(miRNA == mir) %>% pull(Target)
  targets <- targets[targets %in% rownames(target_expr_by_ct)]
  
  if (length(targets) > 0) {
    # 【v5修复】target_expr_by_ct的列名可能与cell_types不一致
    target_cols <- colnames(target_expr_by_ct)
    for (ct in cell_types) {
      actual_ct <- ct
      if (!(ct %in% target_cols)) {
        idx <- which(gsub("-", "_", target_cols) == ct)
        if (length(idx) > 0) actual_ct <- target_cols[idx[1]]
        else next  # 找不到就跳过
      }
      target_signal <- tryCatch(
        mean(target_expr_by_ct[targets, actual_ct, drop=FALSE], na.rm = TRUE),
        error = function(e) { NA }
      )
      if (!is.na(target_signal)) mirna_scores[mir, ct] <- target_signal
    }
  }
}

# --- 综合评分 ---
# 标准化两个矩阵到[0,1]
normalize_matrix <- function(mat) {
  if (all(mat == 0, na.rm = TRUE)) {
    cat("  WARNING: All-zero matrix detected. Returning uniform scores (0.5).\n")
    mat[] <- 0.5
    return(mat)
  }
  mat_range <- range(mat, na.rm = TRUE)
  if (is.na(mat_range[1]) || is.na(mat_range[2]) || mat_range[2] == mat_range[1]) {
    mat[] <- 0.5
    return(mat)
  }
  (mat - mat_range[1]) / (mat_range[2] - mat_range[1])
}

protein_scores_norm <- normalize_matrix(protein_scores)
mirna_scores_norm <- normalize_matrix(mirna_scores)

cat(paste0("  protein_scores_norm range: [", min(protein_scores_norm), ", ", max(protein_scores_norm), "]\n"))
cat(paste0("  mirna_scores_norm range: [", min(mirna_scores_norm), ", ", max(mirna_scores_norm), "]\n"))

# 保存评分矩阵
write.csv(protein_scores_norm,
          file.path(res_dir, "tables", "protein_cargo_celltype_scores.csv"))
write.csv(mirna_scores_norm,
          file.path(res_dir, "tables", "miRNA_cargo_celltype_scores.csv"))

# 2.4 可视化: Cargo × CellType 匹配热图（核心图之一）
pdf(file.path(res_dir, "figures", "04_docking", "cargo_celltype_matching_protein.pdf"),
    width = 12, height = 14)

# 只展示有信号的cargo
active_proteins <- rownames(protein_scores_norm)[rowSums(protein_scores_norm) > 0]
if (length(active_proteins) > 0) {
  Heatmap(
    protein_scores_norm[active_proteins, , drop = FALSE],
    name = "Matching\nScore",
    col = colorRamp2(c(0, 0.5, 1), c("white", "#FFD700", "#E64B35")),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 45,
    row_names_gp = gpar(fontsize = 8),
    column_title = "Protein Cargo → Cell Type Matching",
    row_title = "Exosomal Protein Cargo",
    row_split = cargo_db$proteins$therapeutic_relevance[match(active_proteins, cargo_db$proteins$gene)],
    row_title_rot = 0
  )
}
dev.off()

pdf(file.path(res_dir, "figures", "04_docking", "cargo_celltype_matching_miRNA.pdf"),
    width = 12, height = 10)

active_mirnas <- rownames(mirna_scores_norm)[rowSums(mirna_scores_norm) > 0]
if (length(active_mirnas) > 0) {
  Heatmap(
    mirna_scores_norm[active_mirnas, , drop = FALSE],
    name = "Target\nExpression",
    col = colorRamp2(c(0, 0.5, 1), c("white", "#7E6148", "#3C5488")),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 45,
    row_names_gp = gpar(fontsize = 8),
    column_title = "miRNA Cargo → Cell Type Target Matching",
    row_title = "Exosomal miRNA Cargo"
  )
}
dev.off()

# ============================================================
# PART 3: CellChat 通讯分析 + 外泌体虚拟信号源
# ============================================================

cat("\n--- Part 3: CellChat Communication Analysis ---\n")

# 3.1 常规CellChat分析（正常 vs 损伤）
cat("  Running CellChat analysis...\n")

cellchat_analysis <- function(seurat_obj, condition_label = "all") {
  # 准备数据
  data_input <- GetAssayData(seurat_obj, layer = "data")
  meta <- seurat_obj@meta.data
  meta$labels <- meta$cell_type
  
  cellchat <- createCellChat(object = data_input, meta = meta, group.by = "labels")
  
  # 设置LR数据库
  CellChatDB.use <- subsetDB(CellChatDB.human, 
    search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
  cellchat@DB <- CellChatDB.use
  
  # 分析流程
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  return(cellchat)
}

# 对全部细胞运行
cellchat_all <- cellchat_analysis(merged, "all")

# 如果有条件信息，分别分析
cellchat_list <- list(All = cellchat_all)

if (!is.null(condition_col)) {
  conditions <- unique(merged@meta.data[[condition_col]])
  for (cond in conditions) {
    cond_cells <- subset(merged, cells = colnames(merged)[merged@meta.data[[condition_col]] == cond])
    if (ncol(cond_cells) > 100) {
      cat(paste0("  Running CellChat for condition: ", cond, "\n"))
      cellchat_list[[cond]] <- tryCatch(
        cellchat_analysis(cond_cells, cond),
        error = function(e) { cat(paste0("    Failed: ", e$message, "\n")); NULL }
      )
    }
  }
}

# 3.2 通讯网络可视化
pdf(file.path(res_dir, "figures", "04_docking", "cellchat_network_overview.pdf"),
    width = 12, height = 10)
netVisual_circle(cellchat_all@net$count, 
                  vertex.weight = table(cellchat_all@idents),
                  weight.scale = TRUE, label.edge = FALSE,
                  title.name = "Cell-Cell Communication Network")
dev.off()

# 通路级别通讯
pdf(file.path(res_dir, "figures", "04_docking", "cellchat_pathway_contribution.pdf"),
    width = 14, height = 8)
netAnalysis_contribution(cellchat_all, signaling = cellchat_all@netP$pathways[1:min(20, length(cellchat_all@netP$pathways))])
dev.off()

# 3.3 识别外泌体可修复的受损通讯通路
cat("  Identifying exosome-repairable communication pathways...\n")

# 提取所有活跃的信号通路
active_pathways <- cellchat_all@netP$pathways
cat(paste0("  Active signaling pathways: ", length(active_pathways), "\n"))

# 外泌体cargo覆盖的通路
therapeutic_pathways_flat <- unlist(config$exosome$therapeutic_pathways)

# 找到外泌体cargo可以增强/修复的通路
exosome_targetable_pathways <- active_pathways[
  sapply(active_pathways, function(p) {
    any(sapply(therapeutic_pathways_flat, function(tp) grepl(tp, p, ignore.case = TRUE)))
  })
]

# 补充：直接检查cargo蛋白是否出现在CellChat信号中
cargo_in_cellchat <- cargo_proteins[cargo_proteins %in% 
  unlist(strsplit(c(CellChatDB$interaction$ligand, CellChatDB$interaction$receptor), "_"))]

cat(paste0("  Cargo proteins in CellChatDB: ", length(cargo_in_cellchat), "\n"))
cat(paste0("  Exosome-targetable pathways: ", length(exosome_targetable_pathways), "\n"))

# 3.4 外泌体治疗性通讯重建图
# 概念：将外泌体视为虚拟信号源，预测其对通讯网络的修复效果

# 对每条外泌体可靶向的通路，分析外泌体cargo可以作为哪些配体
exo_repair_map <- data.frame(
  Pathway = character(),
  Cargo_Ligand = character(),
  Target_Receptor = character(),
  Sender_CellType = character(),
  Receiver_CellType = character(),
  Original_Strength = numeric(),
  stringsAsFactors = FALSE
)

lr_pairs <- tryCatch(cellchat_all@LR$LRsig, error = function(e) { 
  cat(paste0("  LR pairs not available: ", e$message, "\n")); NULL 
})

if (!is.null(lr_pairs) && nrow(lr_pairs) > 0) {
  # 过滤掉receptor或ligand为空的行
  valid_lr_idx <- which(
    !is.na(lr_pairs$ligand) & !is.na(lr_pairs$receptor) &
    nchar(as.character(lr_pairs$ligand)) > 0 & nchar(as.character(lr_pairs$receptor)) > 0
  )
  
  for (i in valid_lr_idx) {
    ligands <- unlist(strsplit(as.character(lr_pairs$ligand[i]), "_"))
    receptors <- unlist(strsplit(as.character(lr_pairs$receptor[i]), "_"))
    
    cargo_match <- intersect(ligands, cargo_proteins)
    if (length(cargo_match) > 0) {
      pathway <- if("pathway_name" %in% colnames(lr_pairs)) as.character(lr_pairs$pathway_name[i]) else "Unknown"
      exo_repair_map <- rbind(exo_repair_map, data.frame(
        Pathway = pathway,
        Cargo_Ligand = paste(cargo_match, collapse=","),
        Target_Receptor = paste(receptors, collapse=","),
        Sender_CellType = "Exosome(Virtual)",
        Receiver_CellType = "All_expressing_receptor",
        Original_Strength = 1,
        stringsAsFactors = FALSE
      ))
    }
  }
} else {
  cat("  Using cargo_ligand_matches as repair map fallback.\n")
  if (nrow(cargo_ligand_matches) > 0 && "pathway_name" %in% colnames(cargo_ligand_matches)) {
    for (i in seq_len(nrow(cargo_ligand_matches))) {
      exo_repair_map <- rbind(exo_repair_map, data.frame(
        Pathway = as.character(cargo_ligand_matches$pathway_name[i]),
        Cargo_Ligand = as.character(cargo_ligand_matches$ligand[i]),
        Target_Receptor = as.character(cargo_ligand_matches$receptor[i]),
        Sender_CellType = "Exosome(Virtual)",
        Receiver_CellType = "All_expressing_receptor",
        Original_Strength = 1,
        stringsAsFactors = FALSE
      ))
    }
  }
}

exo_repair_map <- distinct(exo_repair_map)
write.csv(exo_repair_map,
          file.path(res_dir, "tables", "exosome_communication_repair_map.csv"),
          row.names = FALSE)

cat(paste0("  Exosome-mediated communication pairs: ", nrow(exo_repair_map), "\n"))

# ============================================================
# PART 4: 优先级评分与综合排名
# ============================================================

cat("\n--- Part 4: Priority Scoring ---\n")

# 4.1 细胞亚群"可治疗性"评分
cat("  Computing cell type treatability scores...\n")

treatability_scores <- data.frame(cell_type = cell_types)

# 【v5修复】获取receptor_expr和protein_scores的实际列名（Seurat可能把_替换为-）
actual_receptor_cols <- colnames(receptor_expr)
actual_protein_cols <- colnames(protein_scores_norm)
# 建立cell_types到实际列名的映射
ct_to_receptor_col <- sapply(cell_types, function(ct) {
  if (ct %in% actual_receptor_cols) return(ct)
  # 尝试下划线/横线互换
  idx <- which(actual_receptor_cols == gsub("-", "_", ct))
  if (length(idx) > 0) return(actual_receptor_cols[idx[1]])
  idx <- which(gsub("-", "_", actual_receptor_cols) == ct)
  if (length(idx) > 0) return(actual_receptor_cols[idx[1]])
  return(NULL)
}, USE.NAMES = FALSE)

ct_to_protein_col <- sapply(cell_types, function(ct) {
  if (ct %in% actual_protein_cols) return(ct)
  idx <- which(actual_protein_cols == gsub("-", "_", ct))
  if (length(idx) > 0) return(actual_protein_cols[idx[1]])
  idx <- which(gsub("-", "_", actual_protein_cols) == ct)
  if (length(idx) > 0) return(actual_protein_cols[idx[1]])
  return(NULL)
}, USE.NAMES = FALSE)

# 因素1: 受体表达丰富度
valid_rec_rows <- expressed_receptors[expressed_receptors %in% rownames(receptor_expr)]
if (length(valid_rec_rows) > 0 && length(ct_to_receptor_col) > 0) {
  valid_ct_rec <- ct_to_receptor_col[!is.na(ct_to_receptor_col) & ct_to_receptor_col %in% actual_receptor_cols]
  valid_ct_rec <- valid_ct_rec[!duplicated(valid_ct_rec)]
  if (length(valid_ct_rec) > 0) {
    treatability_scores$receptor_richness <- NA
    for (ci in seq_along(cell_types)) {
      if (!is.na(ct_to_receptor_col[ci]) && ct_to_receptor_col[ci] %in% actual_receptor_cols) {
        tryCatch({
          treatability_scores$receptor_richness[ci] <- mean(receptor_expr[valid_rec_rows, ct_to_receptor_col[ci], drop=FALSE], na.rm=TRUE)
        }, error = function(e) { NA })
      }
    }
  } else {
    treatability_scores$receptor_richness <- 0
  }
} else {
  treatability_scores$receptor_richness <- 0
}

# 因素2: 蛋白cargo匹配度
if (length(ct_to_protein_col) > 0) {
  treatability_scores$protein_cargo_match <- NA
  for (ci in seq_along(cell_types)) {
    if (!is.na(ct_to_protein_col[ci]) && ct_to_protein_col[ci] %in% actual_protein_cols) {
      tryCatch({
        treatability_scores$protein_cargo_match[ci] <- mean(protein_scores_norm[, ct_to_protein_col[ci], drop=FALSE], na.rm=TRUE)
      }, error = function(e) { NA })
    }
  }
} else {
  treatability_scores$protein_cargo_match <- 0
}

# 因素3: miRNA靶基因表达
actual_mirna_cols <- colnames(mirna_scores_norm)
ct_to_mirna_col <- sapply(cell_types, function(ct) {
  if (ct %in% actual_mirna_cols) return(ct)
  idx <- which(actual_mirna_cols == gsub("-", "_", ct))
  if (length(idx) > 0) return(actual_mirna_cols[idx[1]])
  idx <- which(gsub("-", "_", actual_mirna_cols) == ct)
  if (length(idx) > 0) return(actual_mirna_cols[idx[1]])
  return(NULL)
}, USE.NAMES = FALSE)

if (length(ct_to_mirna_col) > 0) {
  treatability_scores$mirna_target_expr <- NA
  for (ci in seq_along(cell_types)) {
    if (!is.na(ct_to_mirna_col[ci]) && ct_to_mirna_col[ci] %in% actual_mirna_cols) {
      tryCatch({
        treatability_scores$mirna_target_expr[ci] <- mean(mirna_scores_norm[, ct_to_mirna_col[ci], drop=FALSE], na.rm=TRUE)
      }, error = function(e) { NA })
    }
  }
} else {
  treatability_scores$mirna_target_expr <- 0
}

# 因素4: 细胞通讯参与度（作为receiver的互作数量）
if (!is.null(cellchat_all@net$count)) {
  comm_counts <- cellchat_all@net$count
  # 【v5修复】cell_types与comm_counts的rownames可能_/-不一致
  comm_rownames <- rownames(comm_counts)
  comm_colnames <- colnames(comm_counts)
  treatability_scores$communication_involvement <- 0
  for (ci in seq_along(cell_types)) {
    ct <- cell_types[ci]
    # 尝试多种匹配方式
    ct_in_comm <- NULL
    if (ct %in% comm_colnames) {
      ct_in_comm <- ct
    } else {
      alt <- gsub("-", "_", ct)
      if (alt %in% comm_colnames) ct_in_comm <- alt
      else {
        idx <- which(gsub("-", "_", comm_colnames) == ct)
        if (length(idx) > 0) ct_in_comm <- comm_colnames[idx[1]]
      }
    }
    if (!is.null(ct_in_comm)) {
      tryCatch({
        treatability_scores$communication_involvement[ci] <- sum(comm_counts[, ct_in_comm, drop=FALSE], na.rm = TRUE)
      }, error = function(e) { NA })
    }
  }
}

# 综合评分
score_cols <- c("receptor_richness", "protein_cargo_match", "mirna_target_expr")
if ("communication_involvement" %in% colnames(treatability_scores)) {
  score_cols <- c(score_cols, "communication_involvement")
}

# 标准化各因素到[0,1]
for (col in score_cols) {
  val <- treatability_scores[[col]]
  rng <- range(val, na.rm = TRUE)
  if (rng[2] > rng[1]) {
    treatability_scores[[col]] <- (val - rng[1]) / (rng[2] - rng[1])
  }
}

treatability_scores$composite_score <- rowMeans(treatability_scores[, score_cols], na.rm = TRUE)
treatability_scores <- treatability_scores %>% arrange(desc(composite_score))

write.csv(treatability_scores,
          file.path(res_dir, "tables", "cell_type_treatability_scores.csv"),
          row.names = FALSE)

cat("  Cell type treatability ranking:\n")
print(treatability_scores[, c("cell_type", "composite_score")])

# 4.2 最优cargo组合推荐（按细胞类型）
cat("  Generating optimal cargo recommendations...\n")

cargo_recommendations <- list()
for (ct in head(treatability_scores$cell_type, 5)) {  # Top 5 most treatable
  
  # 【v5修复】映射到实际列名
  actual_ct_prot <- NULL
  if (ct %in% colnames(protein_scores_norm)) actual_ct_prot <- ct
  else {
    idx <- which(gsub("-", "_", colnames(protein_scores_norm)) == ct)
    if (length(idx) > 0) actual_ct_prot <- colnames(protein_scores_norm)[idx[1]]
  }
  
  actual_ct_mirna <- NULL
  if (ct %in% colnames(mirna_scores_norm)) actual_ct_mirna <- ct
  else {
    idx <- which(gsub("-", "_", colnames(mirna_scores_norm)) == ct)
    if (length(idx) > 0) actual_ct_mirna <- colnames(mirna_scores_norm)[idx[1]]
  }
  
  # Best protein cargo for this cell type
  top_proteins <- character(0)
  if (!is.null(actual_ct_prot)) {
    tryCatch({
      top_vec <- sort(protein_scores_norm[, actual_ct_prot, drop=FALSE], decreasing = TRUE)
      top_proteins <- names(top_vec[top_vec > 0])[1:min(10, sum(top_vec > 0))]
    }, error = function(e) { character(0) })
  }
  
  # Best miRNA cargo
  top_mirnas <- character(0)
  if (!is.null(actual_ct_mirna)) {
    tryCatch({
      mirna_vec <- sort(mirna_scores_norm[, actual_ct_mirna, drop=FALSE], decreasing = TRUE)
      top_mirnas <- names(mirna_vec[mirna_vec > 0])[1:min(5, sum(mirna_vec > 0))]
    }, error = function(e) { character(0) })
  }
  
  cargo_recommendations[[ct]] <- list(
    cell_type = ct,
    recommended_proteins = top_proteins,
    recommended_mirnas = top_mirnas,
    treatability_score = treatability_scores$composite_score[treatability_scores$cell_type == ct]
  )
}

# 保存推荐
rec_df <- do.call(rbind, lapply(cargo_recommendations, function(r) {
  data.frame(
    Cell_Type = r$cell_type,
    Treatability_Score = round(r$treatability_score, 3),
    Top_Protein_Cargo = paste(r$recommended_proteins, collapse = "; "),
    Top_miRNA_Cargo = paste(r$recommended_mirnas, collapse = "; "),
    stringsAsFactors = FALSE
  )
}))

write.csv(rec_df,
          file.path(res_dir, "tables", "optimal_cargo_recommendations.csv"),
          row.names = FALSE)

# 4.3 Sankey图: Cargo → Pathway → Cell Type
cat("  Generating Sankey diagram...\n")

# 构建Sankey数据
sankey_data <- data.frame(
  Cargo = character(),
  Pathway = character(), 
  CellType = character(),
  Score = numeric(),
  stringsAsFactors = FALSE
)

for (ct in names(cargo_recommendations)) {
  rec <- cargo_recommendations[[ct]]
  for (prot in rec$recommended_proteins[1:min(5, length(rec$recommended_proteins))]) {
    # 查找该cargo蛋白参与的通路
    cargo_pathways <- exo_repair_map %>%
      filter(grepl(prot, Cargo_Ligand)) %>%
      pull(Pathway) %>%
      unique()
    
    if (length(cargo_pathways) == 0) {
      cargo_pathways <- cargo_db$proteins$therapeutic_relevance[cargo_db$proteins$gene == prot]
    }
    
    for (pw in cargo_pathways[1:min(2, length(cargo_pathways))]) {
      # 【v5修复】安全取值
      score_val <- tryCatch({
        sc <- ct
        if (!(ct %in% colnames(protein_scores_norm))) {
          idx <- which(gsub("-", "_", colnames(protein_scores_norm)) == ct)
          if (length(idx) > 0) sc <- colnames(protein_scores_norm)[idx[1]] else { 0 }
        } else {
          as.numeric(protein_scores_norm[prot, sc])
        }
      }, error = function(e) { 0 })
      sankey_data <- rbind(sankey_data, data.frame(
        Cargo = prot, Pathway = pw, CellType = ct,
        Score = score_val
      ))
    }
  }
}

if (nrow(sankey_data) > 0) {
  sankey_data <- sankey_data %>% filter(Score > 0) %>% distinct()
  
  p_sankey <- ggplot(sankey_data,
    aes(axis1 = Cargo, axis2 = Pathway, axis3 = CellType, y = Score)) +
    geom_alluvium(aes(fill = CellType), alpha = 0.7) +
    geom_stratum(width = 0.3, fill = "grey90", color = "grey50") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.5) +
    scale_fill_npg() +
    theme_minimal() +
    labs(title = "Exosome Cargo → Pathway → Target Cell Type Flow",
         y = "Matching Score") +
    theme(legend.position = "right",
          axis.text.x = element_blank())
  
  pdf(file.path(res_dir, "figures", "04_docking", "sankey_cargo_pathway_celltype.pdf"),
      width = 16, height = 12)
  print(p_sankey)
  dev.off()
}

# ============================================================
# PART 5: 保存所有结果
# ============================================================

cat("\n--- Saving all results ---\n")

docking_results <- list(
  receptor_expr = receptor_expr,
  receptor_zscore = receptor_zscore,
  protein_scores = protein_scores_norm,
  mirna_scores = mirna_scores_norm,
  cellchat = cellchat_list,
  exo_repair_map = exo_repair_map,
  treatability_scores = treatability_scores,
  cargo_recommendations = cargo_recommendations
)

saveRDS(docking_results, file.path(data_dir, "processed", "docking_results.rds"))

# 综合结果Excel
wb <- createWorkbook()
addWorksheet(wb, "Treatability_Scores")
writeData(wb, "Treatability_Scores", treatability_scores)
addWorksheet(wb, "Cargo_Recommendations")
writeData(wb, "Cargo_Recommendations", rec_df)
addWorksheet(wb, "Communication_Repair_Map")
writeData(wb, "Communication_Repair_Map", exo_repair_map)
saveWorkbook(wb, file.path(res_dir, "tables", "molecular_docking_summary.xlsx"),
             overwrite = TRUE)

cat("\n=== Module 3 Complete ===\n")
cat("  All docking results saved.\n")
cat("  Next step: Rscript scripts/04_validation/step1_cross_validation.R\n")
