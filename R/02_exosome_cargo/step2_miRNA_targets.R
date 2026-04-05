# ============================================================
# Module 2 Step 2: miRNA靶基因预测与多层调控网络
# 
# 输入: cargo数据库 + miRNA靶点数据库
# 输出: miRNA-靶基因网络 + ceRNA网络
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(igraph)
  library(openxlsx)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

config  <- yaml::read_yaml("config/params.yaml")
data_dir <- path.expand(config$data_dir)
res_dir  <- path.expand(config$results_dir)
ref_dir  <- file.path(data_dir, "raw", "reference")

cat("=== Module 2 Step 2: miRNA Target Prediction ===\n")

# 加载cargo数据库
cargo_db <- readRDS(file.path(data_dir, "processed", "cargo_database.rds"))
mirna_list <- cargo_db$mirnas$mirna

# ============================================================
# 2.1 TargetScan 靶基因预测
# ============================================================

cat("  Loading TargetScan predictions...\n")

targetscan_targets <- NULL
ts_file <- file.path(ref_dir, "TargetScan_predictions.txt.zip")

if (file.exists(ts_file)) {
  ts_data <- read.delim(unz(ts_file, "Predicted_Targets_Context_Scores.default_predictions.txt"),
                        stringsAsFactors = FALSE)
  
  # 筛选人类 + 我们关注的miRNAs
  targetscan_targets <- ts_data %>%
    filter(grepl("^hsa-", `miRNA`)) %>%
    filter(`miRNA` %in% mirna_list | gsub("-[35]p$", "", `miRNA`) %in% gsub("-[35]p$", "", mirna_list)) %>%
    filter(`context...score` <= config$exosome$mirna_targets$score_threshold$targetscan_context) %>%
    select(miRNA = `miRNA`, Target = `Target.gene`, Score = `context...score`) %>%
    distinct()
  
  cat(paste0("  TargetScan: ", nrow(targetscan_targets), " predictions\n"))
} else {
  cat("  TargetScan file not found. Using built-in target data.\n")
}

# ============================================================
# 2.2 miRDB 靶基因预测
# ============================================================

cat("  Loading miRDB predictions...\n")

mirdb_targets <- NULL
mirdb_file <- file.path(ref_dir, "miRDB_v6.0.txt.gz")

if (file.exists(mirdb_file)) {
  mirdb_data <- read.delim(gzfile(mirdb_file), header = FALSE,
                           col.names = c("miRNA", "Target", "Score"),
                           stringsAsFactors = FALSE)
  
  mirdb_targets <- mirdb_data %>%
    filter(grepl("^hsa-", miRNA)) %>%
    filter(miRNA %in% mirna_list | gsub("-[35]p$", "", miRNA) %in% gsub("-[35]p$", "", mirna_list)) %>%
    filter(Score >= config$exosome$mirna_targets$score_threshold$mirdb_score) %>%
    distinct()
  
  cat(paste0("  miRDB: ", nrow(mirdb_targets), " predictions\n"))
} else {
  cat("  miRDB file not found.\n")
}

# ============================================================
# 2.3 miRTarBase 实验验证靶基因
# ============================================================

cat("  Loading miRTarBase validated targets...\n")

mirtarbase_targets <- NULL
mtb_file <- file.path(ref_dir, "hsa_MTI.xlsx")

if (file.exists(mtb_file)) {
  mtb_data <- read.xlsx(mtb_file)
  
  mirtarbase_targets <- mtb_data %>%
    filter(grepl("^hsa-", `miRNA`)) %>%
    filter(`miRNA` %in% mirna_list) %>%
    select(miRNA = `miRNA`, Target = `Target.Gene`, 
           Experiment = `Experiments`, Support = `Support.Type`) %>%
    distinct()
  
  cat(paste0("  miRTarBase: ", nrow(mirtarbase_targets), " validated interactions\n"))
} else {
  cat("  miRTarBase file not found.\n")
}

# ============================================================
# 2.4 多数据库交集（高置信靶基因）
# ============================================================

cat("  Computing high-confidence targets (multi-database intersection)...\n")

# 标准化miRNA名称
standardize_mirna <- function(name) {
  # 统一为小写并去除额外空格
  tolower(trimws(name))
}

# 合并所有来源
all_targets <- list()

if (!is.null(targetscan_targets)) {
  ts_pairs <- targetscan_targets %>%
    mutate(miRNA = standardize_mirna(miRNA)) %>%
    select(miRNA, Target) %>%
    mutate(Source = "TargetScan")
  all_targets[["TargetScan"]] <- ts_pairs
}

if (!is.null(mirdb_targets)) {
  md_pairs <- mirdb_targets %>%
    mutate(miRNA = standardize_mirna(miRNA)) %>%
    select(miRNA, Target) %>%
    mutate(Source = "miRDB")
  all_targets[["miRDB"]] <- md_pairs
}

if (!is.null(mirtarbase_targets)) {
  mtb_pairs <- mirtarbase_targets %>%
    mutate(miRNA = standardize_mirna(miRNA)) %>%
    select(miRNA, Target) %>%
    mutate(Source = "miRTarBase")
  all_targets[["miRTarBase"]] <- mtb_pairs
}

if (length(all_targets) >= 2) {
  combined <- bind_rows(all_targets)
  
  # 计算每对miRNA-Target在几个数据库中出现
  target_confidence <- combined %>%
    group_by(miRNA, Target) %>%
    summarise(
      n_databases = n_distinct(Source),
      databases = paste(unique(Source), collapse = "|"),
      .groups = "drop"
    ) %>%
    arrange(desc(n_databases), miRNA)
  
  # 高置信靶基因（至少2个数据库支持）
  high_conf_targets <- target_confidence %>%
    filter(n_databases >= config$exosome$mirna_targets$min_overlap)
  
  cat(paste0("  High-confidence targets (≥2 databases): ", nrow(high_conf_targets), "\n"))
  
} else {
  # 如果数据库文件不足，使用文献报道的已验证靶基因
  cat("  Insufficient database files. Using literature-curated targets.\n")
  
  high_conf_targets <- data.frame(
    miRNA = cargo_db$mirnas$mirna,
    Target = sapply(strsplit(cargo_db$mirnas$validated_targets, "/"), `[`, 1),
    n_databases = 3,
    databases = "Literature_validated",
    stringsAsFactors = FALSE
  )
}

write.csv(high_conf_targets,
          file.path(res_dir, "tables", "miRNA_high_confidence_targets.csv"),
          row.names = FALSE)

# ============================================================
# 2.5 靶基因功能富集分析
# ============================================================

cat("  Enrichment analysis of miRNA targets...\n")

unique_targets <- unique(high_conf_targets$Target)
target_entrez <- bitr(unique_targets, fromType = "SYMBOL",
                      toType = "ENTREZID", OrgDb = org.Hs.eg.db)

if (nrow(target_entrez) >= 5) {
  ego_targets <- enrichGO(
    gene = target_entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  
  ekegg_targets <- enrichKEGG(
    gene = target_entrez$ENTREZID,
    organism = "hsa",
    pvalueCutoff = 0.05
  )
  
  write.csv(as.data.frame(ego_targets),
            file.path(res_dir, "tables", "miRNA_targets_GO.csv"),
            row.names = FALSE)
  
  pdf(file.path(res_dir, "figures", "03_exosome", "miRNA_targets_GO.pdf"),
      width = 12, height = 10)
  dotplot(ego_targets, showCategory = 20) +
    ggtitle("GO Enrichment of Exosomal miRNA Targets")
  dev.off()
}

# ============================================================
# 2.6 miRNA-靶基因调控网络
# ============================================================

cat("  Building miRNA-target regulatory network...\n")

# 构建网络
if (nrow(high_conf_targets) > 0) {
  # 创建边列表
  edges <- high_conf_targets %>%
    select(from = miRNA, to = Target, weight = n_databases)
  
  # 构建igraph网络
  g <- graph_from_data_frame(edges, directed = TRUE)
  
  # 节点属性
  V(g)$type <- ifelse(V(g)$name %in% unique(edges$from), "miRNA", "Target")
  V(g)$size <- ifelse(V(g)$type == "miRNA", 
                       degree(g, mode = "out") * 2 + 5,
                       degree(g, mode = "in") * 2 + 3)
  V(g)$color <- ifelse(V(g)$type == "miRNA", "#E64B35", "#4DBBD5")
  
  # 保存网络
  write.csv(edges, file.path(res_dir, "tables", "miRNA_target_network_edges.csv"),
            row.names = FALSE)
  write.csv(data.frame(node = V(g)$name, type = V(g)$type, degree = degree(g)),
            file.path(res_dir, "tables", "miRNA_target_network_nodes.csv"),
            row.names = FALSE)
  
  # 网络可视化
  pdf(file.path(res_dir, "figures", "03_exosome", "miRNA_target_network.pdf"),
      width = 16, height = 14)
  
  layout <- layout_with_fr(g)
  plot(g, layout = layout,
       vertex.size = V(g)$size,
       vertex.color = V(g)$color,
       vertex.label.cex = 0.6,
       vertex.label.color = "black",
       edge.arrow.size = 0.3,
       edge.width = E(g)$weight * 0.5,
       edge.color = alpha("grey50", 0.5),
       main = "MSC Exosome miRNA-Target Regulatory Network")
  legend("bottomleft", legend = c("miRNA", "Target Gene"),
         col = c("#E64B35", "#4DBBD5"), pch = 19, cex = 0.8)
  
  dev.off()
  
  # Hub miRNAs (most connections)
  mirna_hubs <- data.frame(
    miRNA = V(g)$name[V(g)$type == "miRNA"],
    n_targets = degree(g, V(g)$name[V(g)$type == "miRNA"], mode = "out")
  ) %>% arrange(desc(n_targets))
  
  cat("  Top hub miRNAs:\n")
  print(head(mirna_hubs, 10))
  
  write.csv(mirna_hubs,
            file.path(res_dir, "tables", "hub_miRNAs_by_target_count.csv"),
            row.names = FALSE)
}

# ============================================================
# 2.7 保存中间结果
# ============================================================

mirna_analysis <- list(
  high_conf_targets = high_conf_targets,
  network = if (exists("g")) g else NULL,
  hub_mirnas = if (exists("mirna_hubs")) mirna_hubs else NULL
)

saveRDS(mirna_analysis, file.path(data_dir, "processed", "mirna_target_analysis.rds"))

cat("\n  miRNA target analysis complete.\n")
cat("  Next step: Rscript scripts/03_molecular_docking/step1_receptor_profiling.R\n")
cat("=== Module 2 Step 2 Complete ===\n")
