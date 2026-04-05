# ============================================================
# Step 2: 损伤微环境通讯网络差异分析 (CellChat)
#
# 核心逻辑: 对比 Normal vs Degenerated 组的细胞通讯网络,
#          识别退变状态下特有的"异常通讯通路",
#          重点关注炎症细胞→祖细胞/纤维软骨细胞的致病通路
#
# 输入: data/processed/meniscus_annotated.rds
# 输出: results/figures/03_cellchat_differential/
#       results/tables/cellchat_*.csv
#
# 运行方式: Rscript step2_cellchat_differential.R
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(yaml)
  library(tidyverse)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(igraph)
  library(ggsci)
  library(ggpubr)
  library(ggrepel)
})

config  <- yaml::read_yaml("config/params.yaml")
data_dir <- path.expand(config$data_dir)
res_dir  <- path.expand(config$results_dir)

fig_dir <- file.path(res_dir, "figures", "03_cellchat_differential")
tbl_dir <- file.path(res_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Step 2: 损伤微环境通讯网络差异分析\n")
cat("========================================\n\n")

# ============================================================
# 1. 加载数据
# ============================================================
cat("[1/9] 加载Seurat对象...\n")
obj <- readRDS(file.path(data_dir, "processed", "meniscus_annotated.rds"))

# 检查condition信息
if (!"condition" %in% colnames(obj@meta.data)) {
  stop("错误: Seurat对象中没有 'condition' 列 (Normal/Degenerated)")
}

conditions <- unique(obj$condition)
cat(paste0("  条件分组: ", paste(conditions, collapse=", "), "\n"))

# 检查每组细胞数
for (cond in conditions) {
  n <- sum(obj$condition == cond)
  cat(paste0("  ", cond, ": ", n, " 个细胞\n"))
  if (n < 100) {
    cat(paste0("  ⚠️ 警告: ", cond, " 组细胞数 < 100, CellChat可能不稳定\n"))
  }
}
cat("\n")

# ============================================================
# 2. CellChat准备: 合并小簇
# ============================================================
cat("[2/9] 准备CellChat输入...\n")

# 对于每种细胞类型, 检查在各条件下的细胞数
cell_counts <- table(obj$condition, obj$cell_type)
cat("  各组各细胞类型细胞数:\n")
print(cell_counts)
cat("\n")

# 合并细胞数太少 (< 10个) 的簇为 "Other"
min_cells_per_type <- 10
type_counts <- table(obj$cell_type)
small_types <- names(type_counts[type_counts < min_cells_per_type])

if (length(small_types) > 0) {
  cat(paste0("  合并小簇 (< ", min_cells_per_type, " 个细胞): ",
             paste(small_types, collapse=", "), " → Other\n"))
  obj$cell_type_merged <- as.character(obj$cell_type)
  obj$cell_type_merged[obj$cell_type %in% small_types] <- "Other"
  group_by_col <- "cell_type_merged"
} else {
  obj$cell_type_merged <- obj$cell_type
  group_by_col <- "cell_type"
}

# ============================================================
# 3. 分别创建 Normal 和 Degenerated 的 CellChat 对象
# ============================================================
cat("[3/9] 创建CellChat对象...\n")

# 通用CellChat分析函数
run_cellchat <- function(seurat_subset, label) {
  cat(paste0("  构建 ", label, " CellChat对象...\n"))
  
  # 准备数据
  data_input <- GetAssayData(seurat_subset, layer = "data")
  meta <- seurat_subset@meta.data
  meta$labels <- meta[[group_by_col]]
  
  cellchat <- createCellChat(object = data_input, meta = meta, group.by = "labels")
  
  # 设置配体-受体数据库
  CellChatDB.use <- subsetDB(CellChatDB.human,
                              search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
  cellchat@DB <- CellChatDB.use
  
  # 标准分析流程
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  cat(paste0("    细胞类型数: ", length(unique(cellchat@idents)), "\n"))
  cat(paste0("    通讯对数: ", nrow(cellchat@net$count), "\n"))
  cat(paste0("    活跃通路数: ", length(cellchat@netP$pathways), "\n"))
  
  return(cellchat)
}

# 分别运行
cellchat_list <- list()
for (cond in conditions) {
  cond_cells <- colnames(obj)[obj$condition == cond]
  obj_cond <- subset(obj, cells = cond_cells)
  
  if (ncol(obj_cond) < 100) {
    cat(paste0("  跳过 ", cond, " (细胞数不足)\n"))
    next
  }
  
  cellchat_list[[cond]] <- tryCatch(
    run_cellchat(obj_cond, cond),
    error = function(e) {
      cat(paste0("  ❌ ", cond, " CellChat失败: ", e$message, "\n"))
      return(NULL)
    }
  )
}

# 移除失败的
cellchat_list <- cellchat_list[!sapply(cellchat_list, is.null)]
cat(paste0("\n  成功创建 ", length(cellchat_list), " 个CellChat对象\n\n"))

# ============================================================
# 4. 通讯数量和强度对比
# ============================================================
cat("[4/9] 对比Normal vs Degenerated通讯差异...\n")

if (length(cellchat_list) >= 2) {
  # 统一细胞类型标签 (确保两组可比)
  cellchat_list <- liftCellChat(cellchat_list)
  
  # 对比通讯数量
  num_comparison <- compareInteractions(cellchat_list, show.legend = TRUE)
  
  pdf(file.path(fig_dir, "Fig_comparison_interaction_number.pdf"), width = 5, height = 5)
  print(num_comparison$count + ggtitle("Number of Interactions"))
  dev.off()
  
  pdf(file.path(fig_dir, "Fig_comparison_interaction_weight.pdf"), width = 5, height = 5)
  print(num_comparison$weight + ggtitle("Interaction Strength (Weight)"))
  dev.off()
  
  # 差异通讯数量和强度
  diff_num <- diffInteraction(cellchat_list, 
                               comparison = c(1, 2),  # 默认: Normal vs Degenerated
                               raw.use = FALSE)
  pdf(file.path(fig_dir, "Fig_diff_interaction_number.pdf"), width = 14, height = 6)
  print(diff_num$count + ggtitle("Differential Interaction Number (Degenerated - Normal)"))
  dev.off()
  
  diff_weight <- diffInteraction(cellchat_list,
                                  comparison = c(1, 2),
                                  raw.use = FALSE)
  pdf(file.path(fig_dir, "Fig_diff_interaction_weight.pdf"), width = 14, height = 6)
  print(diff_weight$weight + ggtitle("Differential Interaction Weight (Degenerated - Normal)"))
  dev.off()
}

# ============================================================
# 5. 信号通路层面差异分析
# ============================================================
cat("[5/9] 信号通路差异分析...\n")

if (length(cellchat_list) >= 2) {
  # 提取所有通路的差异信息
  pathway_diff <- list()
  
  all_pathways <- unique(unlist(lapply(cellchat_list, function(cc) cc@netP$pathways)))
  cat(paste0("  总信号通路: ", length(all_pathways), "\n"))
  
  for (pw in all_pathways) {
    pw_info <- data.frame(Pathway = pw)
    
    for (cond_name in names(cellchat_list)) {
      cc <- cellchat_list[[cond_name]]
      
      if (pw %in% cc@netP$pathways) {
        df <- cc@netP$prob[pw]
        pw_info[[paste0(cond_name, "_prob")]] <- sum(df)
        pw_info[[paste0(cond_name, "_pair_count")]] <- nrow(df)
        
        # 通路级别信息
        pw_df <- cc@netP$pathways[pw]
        pw_info[[paste0(cond_name, "_pvalue")]] <- ifelse(
          "pvalue" %in% colnames(pw_df), pw_df$pvalue, NA
        )
      } else {
        pw_info[[paste0(cond_name, "_prob")]] <- 0
        pw_info[[paste0(cond_name, "_pair_count")]] <- 0
        pw_info[[paste0(cond_name, "_pvalue")]] <- NA
      }
    }
    
    pathway_diff[[pw]] <- pw_info
  }
  
  pathway_diff_df <- bind_rows(pathway_diff)
  
  # 计算变化倍数 (Degenerated / Normal)
  cond_names <- names(cellchat_list)
  if (length(cond_names) >= 2) {
    normal_prob_col <- paste0(cond_names[1], "_prob")
    degen_prob_col  <- paste0(cond_names[2], "_prob")
    
    pathway_diff_df$Fold_Change <- ifelse(
      pathway_diff_df[[normal_prob_col]] > 0,
      pathway_diff_df[[degen_prob_col]] / pathway_diff_df[[normal_prob_col]],
      Inf
    )
    pathway_diff_df$Diff_Prob <- pathway_diff_df[[degen_prob_col]] - pathway_diff_df[[normal_prob_col]]
    pathway_diff_df$Direction <- ifelse(
      pathway_diff_df$Diff_Prob > 0, "Enhanced_in_Degenerated", "Reduced_in_Degenerated"
    )
  }
  
  # 排序并保存
  pathway_diff_df <- pathway_diff_df %>%
    arrange(desc(abs(Diff_Prob)))
  
  write.csv(pathway_diff_df,
            file.path(tbl_dir, "cellchat_pathway_differential.csv"),
            row.names = FALSE)
  
  # 显著增强的通路 (退变相关)
  enhanced_pws <- pathway_diff_df %>%
    filter(Direction == "Enhanced_in_Degenerated", Diff_Prob > 0) %>%
    arrange(desc(Diff_Prob))
  
  cat(paste0("  退变增强的通路: ", nrow(enhanced_pws), "\n"))
  if (nrow(enhanced_pws) > 0) {
    for (i in 1:min(10, nrow(enhanced_pws))) {
      cat(paste0("    ", i, ". ", enhanced_pws$Pathway[i],
                 " (FC=", round(enhanced_pws$Fold_Change[i], 2),
                 ", Δ=", round(enhanced_pws$Diff_Prob[i], 4), ")\n"))
    }
  }
  
  # 显著减弱的通路 (正常保护性)
  reduced_pws <- pathway_diff_df %>%
    filter(Direction == "Reduced_in_Degenerated", Diff_Prob < 0) %>%
    arrange(Diff_Prob)
  
  cat(paste0("  退变减弱的通路: ", nrow(reduced_pws), "\n\n"))
}

# ============================================================
# 6. 关键致病通路的深入分析
# ============================================================
cat("[6/9] 关键致病通路深入分析...\n")

if (length(cellchat_list) >= 2) {
  # 定义关注的关键通路
  focus_pathways <- c(
    "TNF", "IL1", "IL6", "IL17",    # 炎症
    "TGFb", "BMP", "WNT",           # 修复/分化
    "FGF", "PDGF", "VEGF",           # 生长因子
    "CXCL", "CCL",                   # 趋化因子
    "COMPLEMENT", "MACROPHAGE",       # 免疫
    "COLLAGEN", "FN1", "LAMININ",     # ECM
    "THBS", "MIF", "GALECTIN",        # 其他
    "ANGIO", "APP", "SEMaphorin"      # 血管/信号
  )
  
  # 找到在数据中实际存在的通路
  available_focus <- focus_pathways[focus_pathways %in% all_pathways]
  cat(paste0("  关注的关键通路 (存在于数据中): ", length(available_focus), "\n"))
  
  # 对每条关键通路, 绘制 Normal vs Degenerated 对比图
  for (pw in available_focus) {
    tryCatch({
      # 检查通路是否在两组中都存在
      pw_in_groups <- sapply(names(cellchat_list), function(cn) pw %in% cellchat_list[[cn]]@netP$pathways)
      
      if (sum(pw_in_groups) == 0) next
      
      # 通路级别的和弦图对比
      pdf(file.path(fig_dir, paste0("Fig_pathway_", pw, "_normal.pdf")), width = 8, height = 7)
      if (pw_in_groups[1]) {
        cc <- cellchat_list[[names(cellchat_list)[1]]]
        if (pw %in% cc@netP$pathways) {
          netVisual_aggregate(cc, signaling = pw, layout = "circle",
                              title = paste(pw, "- Normal"),
                              vertex.weight = table(cc@idents),
                              weight.scale = TRUE)
        }
      }
      dev.off()
      
      pdf(file.path(fig_dir, paste0("Fig_pathway_", pw, "_degenerated.pdf")), width = 8, height = 7)
      if (pw_in_groups[2]) {
        cc <- cellchat_list[[names(cellchat_list)[2]]]
        if (pw %in% cc@netP$pathways) {
          netVisual_aggregate(cc, signaling = pw, layout = "circle",
                              title = paste(pw, "- Degenerated"),
                              vertex.weight = table(cc@idents),
                              weight.scale = TRUE)
        }
      }
      dev.off()
      
      # 热图: 细胞对之间的通讯概率
      pdf(file.path(fig_dir, paste0("Fig_pathway_", pw, "_heatmap_comparison.pdf")), width = 16, height = 7)
      tryCatch({
        p <- netVisual_heatmap(cellchat_list, signaling = pw,
                               measure = "weight", color.heatmap = "Reds",
                               title = paste(pw, "Communication Probability"))
        print(p)
      }, error = function(e) {
        # 如果对比失败, 单独画
        for (cn in names(cellchat_list)) {
          cc <- cellchat_list[[cn]]
          if (pw %in% cc@netP$pathways) {
            netVisual_heatmap(cc, signaling = pw, measure = "weight",
                              title = paste(pw, "-", cn), color.heatmap = "Reds")
          }
        }
      })
      dev.off()
      
    }, error = function(e) {
      cat(paste0("    通路 ", pw, " 可视化失败: ", e$message, "\n"))
    })
  }
  
  cat("\n")
}

# ============================================================
# 7. 退变特异的炎症→祖细胞/软骨细胞通讯网络
# ============================================================
cat("[7/9] 分析炎症→靶细胞的致病通讯...\n")

if (length(cellchat_list) >= 2) {
  # 定义炎症源细胞和靶细胞
  inflammatory_senders <- c("Macrophage", "Male_Inflammatory", "FC-Inflammatory")
  target_receivers <- c("Meniscus_Progenitor", "Stem_Progenitor",
                        "FC-ECM", "FC-Regulatory", "FC-Degenerated",
                        "Hypertrophic_Chondrocyte", "Secretory_Chondrocyte")
  
  available_senders <- intersect(inflammatory_senders, unique(obj[[group_by_col]]))
  available_receivers <- intersect(target_receivers, unique(obj[[group_by_col]]))
  
  cat(paste0("  炎症源细胞: ", paste(available_senders, collapse=", "), "\n"))
  cat(paste0("  靶细胞: ", paste(available_receivers, collapse=", "), "\n"))
  
  # 提取炎症→靶细胞的通讯
  sender_receiver_comm <- list()
  
  for (cond_name in names(cellchat_list)) {
    cc <- cellchat_list[[cond_name]]
    lr_df <- cc@LR$LRsig
    
    for (sender in available_senders) {
      for (receiver in available_receivers) {
        # CellChat的source/target是group number
        sender_id <- which(levels(cc@idents) == sender)
        receiver_id <- which(levels(cc@idents) == receiver)
        
        if (length(sender_id) == 0 || length(receiver_id) == 0) next
        
        pw_comm <- subsetCommunication(cc,
                                        source.use = sender_id,
                                        target.use = receiver_id)
        
        if (nrow(pw_comm) > 0) {
          pw_comm$Condition <- cond_name
          pw_comm$Sender <- sender
          pw_comm$Receiver <- receiver
          sender_receiver_comm[[paste(cond_name, sender, receiver, sep = "_")]] <- pw_comm
        }
      }
    }
  }
  
  if (length(sender_receiver_comm) > 0) {
    sr_df <- bind_rows(sender_receiver_comm)
    
    # 按通路汇总
    sr_summary <- sr_df %>%
      group_by(Condition, Sender, Receiver, pathway_name) %>%
      summarise(
        Total_prob = sum(prob),
        Interaction_count = n(),
        .groups = "drop"
      )
    
    write.csv(sr_summary,
              file.path(tbl_dir, "cellchat_inflammation_target_communication.csv"),
              row.names = FALSE)
    
    # 找到退变中显著增强的炎症→靶细胞通路
    if (length(cond_names) >= 2) {
      enhanced_inflam <- sr_summary %>%
        pivot_wider(names_from = Condition, values_from = c(Total_prob, Interaction_count),
                    names_sep = "_") %>%
        mutate(
          Diff = !!sym(paste0("Total_prob_", cond_names[2])) - !!sym(paste0("Total_prob_", cond_names[1])),
          FC = ifelse(!!sym(paste0("Total_prob_", cond_names[1])) > 0,
                      !!sym(paste0("Total_prob_", cond_names[2])) / !!sym(paste0("Total_prob_", cond_names[1])),
                      Inf)
        ) %>%
        filter(Diff > 0) %>%
        arrange(desc(Diff))
      
      write.csv(enhanced_inflam,
                file.path(tbl_dir, "cellchat_enhanced_inflammation_pathways.csv"),
                row.names = FALSE)
      
      cat("\n  退变中显著增强的炎症→靶细胞通讯:\n")
      if (nrow(enhanced_inflam) > 0) {
        for (i in 1:min(15, nrow(enhanced_inflam))) {
          cat(paste0("    ", enhanced_inflam$Sender[i], " → ",
                     enhanced_inflam$Receiver[i], " [",
                     enhanced_inflam$pathway_name[i], "] ",
                     "Δprob=", round(enhanced_inflam$Diff[i], 4),
                     " (FC=", round(enhanced_inflam$FC[i], 2), ")\n"))
        }
      } else {
        cat("    未发现显著增强的炎症→靶细胞通讯\n")
      }
    }
  }
  cat("\n")
}

# ============================================================
# 8. 退变"异常通讯网络"综合可视化
# ============================================================
cat("[8/9] 绘制退变异常通讯网络综合图...\n")

if (length(cellchat_list) >= 2) {
  # --- 8.1 合并网络的弦图 (核心图) ---
  pdf(file.path(fig_dir, "Fig_network_overview_normal.pdf"), width = 8, height = 7)
  cc_normal <- cellchat_list[[cond_names[1]]]
  netVisual_circle(cc_normal@net$count,
                   vertex.weight = table(cc_normal@idents),
                   weight.scale = TRUE,
                   label.edge = FALSE,
                   title.name = paste("Normal - Communication Network"))
  dev.off()
  
  pdf(file.path(fig_dir, "Fig_network_overview_degenerated.pdf"), width = 8, height = 7)
  cc_degen <- cellchat_list[[cond_names[2]]]
  netVisual_circle(cc_degen@net$count,
                   vertex.weight = table(cc_degen@idents),
                   weight.scale = TRUE,
                   label.edge = FALSE,
                   title.name = paste("Degenerated - Communication Network"))
  dev.off()
  
  # --- 8.2 差异通讯网络 (NetView) ---
  pdf(file.path(fig_dir, "Fig_diff_network_netview.pdf"), width = 16, height = 7)
  tryCatch({
    # 只取退变中增强的通讯对
    pdf(file.path(fig_dir, "Fig_diff_network_enhanced.pdf"), width = 16, height = 7)
    netVisual_aggregate(cellchat_list,
                        signaling = enhanced_pws$Pathway[1:min(5, nrow(enhanced_pws))],
                        layout = "circle",
                        vertex.weight = table(cellchat_list[[cond_names[2]]]@idents),
                        weight.scale = TRUE,
                        title = "Degeneration-Enhanced Signaling")
    dev.off()
  }, error = function(e) {
    cat(paste0("  NetView可视化跳过: ", e$message, "\n"))
  })
  
  # --- 8.3 信号通路贡献比例对比 ---
  pdf(file.path(fig_dir, "Fig_pathway_contribution_comparison.pdf"), width = 18, height = 8)
  tryCatch({
    pathways_to_show <- intersect(
      cellchat_list[[cond_names[1]]]@netP$pathways,
      cellchat_list[[cond_names[2]]]@netP$pathways
    )
    pathways_to_show <- head(pathways_to_show, 15)
    
    if (length(pathways_to_show) > 0) {
      netAnalysis_contribution(cellchat_list, signaling = pathways_to_show)
    }
  }, error = function(e) {
    cat(paste0("  通路贡献对比跳过: ", e$message, "\n"))
  })
  dev.off()
  
  # --- 8.4 关键信号通路的Gene Expression散点图 ---
  # TNF信号通路为例
  key_signaling_genes <- list(
    TNF    = c("TNF", "TNFRSF1A", "TNFRSF1B"),
    IL1    = c("IL1A", "IL1B", "IL1R1"),
    TGFb   = c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2"),
    BMP    = c("BMP2", "BMP4", "BMP7", "BMPR1A", "BMPR2"),
    IL6    = c("IL6", "IL6R", "IL6ST"),
    VEGF   = c("VEGFA", "VEGFB", "KDR", "FLT1"),
    WNT    = c("WNT5A", "WNT3A", "FZD1", "FZD2", "LRP5")
  )
  
  # 找到在退变中显著增强的通路, 检查其配体/受体表达
  if (nrow(enhanced_pws) > 0) {
    top_enhanced_pws <- head(enhanced_pws$Pathway, 5)
    
    # 在CellChatDB中查找这些通路的配体/受体
    for (pw_name in top_enhanced_pws) {
      matching_pathways <- CellChatDB$interaction$pathway_name[
        grepl(pw_name, CellChatDB$interaction$pathway_name, ignore.case = TRUE)
      ]
      
      if (length(matching_pathways) > 0) {
        pw_genes <- unique(c(
          unlist(strsplit(CellChatDB$interaction$ligand[
            CellChatDB$interaction$pathway_name == matching_pathways[1]], "_")),
          unlist(strsplit(CellChatDB$interaction$receptor[
            CellChatDB$interaction$pathway_name == matching_pathways[1]], "_"))
        ))
        pw_genes <- pw_genes[pw_genes %in% rownames(obj)]
        
        if (length(pw_genes) > 0) {
          pdf(file.path(fig_dir, paste0("Fig_gene_expr_", gsub("[^a-zA-Z0-9]", "_", pw_name), ".pdf")),
              width = 14, height = 6)
          tryCatch({
            p1 <- VlnPlot(obj, features = head(pw_genes, 6), group.by = "condition",
                          pt.size = 0, cols = c("#4DBBD5", "#E64B35")) +
              RotatedAxis() + ggtitle(paste(pw_name, "- Gene Expression: Normal vs Degenerated"))
            print(p1)
          }, error = function(e) NULL)
          dev.off()
        }
      }
    }
  }
  
  # --- 8.5 综合热图: 通路 × 细胞对 通讯概率对比 ---
  pdf(file.path(fig_dir, "Fig_pathway_cellpair_heatmap_comparison.pdf"), width = 16, height = 10)
  tryCatch({
    netVisual_heatmap(cellchat_list,
                      measure = "weight",
                      color.heatmap = "Reds",
                      title = "Signaling Pathway Comparison: Normal vs Degenerated")
  }, error = function(e) {
    cat(paste0("  综合热图跳过: ", e$message, "\n"))
  })
  dev.off()
}

# ============================================================
# 9. 保存结果与摘要
# ============================================================
cat("[9/9] 保存结果...\n")

# 保存CellChat对象
saveRDS(cellchat_list, file.path(data_dir, "processed", "cellchat_differential_list.rds"))

# 生成综合报告
report <- data.frame(
  Item = character(),
  Value = character(),
  stringsAsFactors = FALSE
)

if (length(cellchat_list) >= 2) {
  for (cn in names(cellchat_list)) {
    cc <- cellchat_list[[cn]]
    report <- rbind(report, data.frame(
      Item = paste0(cn, "_cell_types"),
      Value = as.character(length(unique(cc@idents))),
      stringsAsFactors = FALSE
    ))
    report <- rbind(report, data.frame(
      Item = paste0(cn, "_active_pathways"),
      Value = as.character(length(cc@netP$pathways)),
      stringsAsFactors = FALSE
    ))
    report <- rbind(report, data.frame(
      Item = paste0(cn, "_interaction_pairs"),
      Value = as.character(nrow(cc@net$count)),
      stringsAsFactors = FALSE
    ))
  }
  
  report <- rbind(report, data.frame(
    Item = "Enhanced_pathways_in_degeneration",
    Value = as.character(nrow(enhanced_pws)),
    stringsAsFactors = FALSE
  ))
  report <- rbind(report, data.frame(
    Item = "Reduced_pathways_in_degeneration",
    Value = as.character(nrow(reduced_pws)),
    stringsAsFactors = FALSE
  ))
  
  report <- rbind(report, data.frame(
    Item = "Top_enhanced_pathway",
    Value = ifelse(nrow(enhanced_pws) > 0, enhanced_pws$Pathway[1], "NA"),
    stringsAsFactors = FALSE
  ))
  report <- rbind(report, data.frame(
    Item = "Top_reduced_pathway",
    Value = ifelse(nrow(reduced_pws) > 0, reduced_pws$Pathway[1], "NA"),
    stringsAsFactors = FALSE
  ))
}

write.csv(report,
          file.path(tbl_dir, "cellchat_differential_summary.csv"),
          row.names = FALSE)

# 打印摘要
cat("\n========================================\n")
cat("分析完成! 核心结果摘要:\n")
cat("========================================\n")

if (length(cellchat_list) >= 2) {
  cat(paste0("  Normal组活跃通路: ", length(cellchat_list[[cond_names[1]]]@netP$pathways), "\n"))
  cat(paste0("  Degenerated组活跃通路: ", length(cellchat_list[[cond_names[2]]]@netP$pathways), "\n"))
  cat(paste0("  退变增强的通路: ", nrow(enhanced_pws), "\n"))
  cat(paste0("  退变减弱的通路: ", nrow(reduced_pws), "\n"))
  
  if (nrow(enhanced_pws) > 0) {
    cat("\n  🔴 退变增强的Top通路:\n")
    for (i in 1:min(10, nrow(enhanced_pws))) {
      cat(paste0("    ", i, ". ", enhanced_pws$Pathway[i],
                 " (FC=", round(enhanced_pws$Fold_Change[i], 2), ")\n"))
    }
  }
  
  if (nrow(reduced_pws) > 0) {
    cat("\n  🔵 退变减弱的Top通路 (正常保护性):\n")
    for (i in 1:min(10, nrow(reduced_pws))) {
      cat(paste0("    ", i, ". ", reduced_pws$Pathway[i],
                 " (FC=", round(reduced_pws$Fold_Change[i], 2), ")\n"))
    }
  }
}

cat("\n输出文件:\n")
cat(paste0("  图表: ", fig_dir, "\n"))
cat(paste0("  表格: ", tbl_dir, "/cellchat_*.csv\n"))
cat(paste0("  CellChat对象: data/processed/cellchat_differential_list.rds\n"))
cat("\n========================================\n")
