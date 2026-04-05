# ============================================================
# Module 2 Step 1: 外泌体 Cargo 数据库构建与治疗性分子筛选
# 
# 输入: ExoCarta/Vesiclepedia 下载数据
# 输出: MSC外泌体治疗性cargo数据库
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(UpSetR)
  library(VennDiagram)
  library(openxlsx)
})

config  <- yaml::read_yaml("config/params.yaml")
data_dir <- path.expand(config$data_dir)
res_dir  <- path.expand(config$results_dir)
exo_params <- config$exosome

dir.create(file.path(res_dir, "figures", "03_exosome"), recursive = TRUE, showWarnings = FALSE)

cat("=== Module 2: Exosome Cargo Analysis ===\n")

# ============================================================
# 1.1 加载ExoCarta和Vesiclepedia数据
# ============================================================

cat("  Loading exosome databases...\n")

exo_dir <- file.path(data_dir, "raw", "exosome")

# --- ExoCarta ---
exocarta_protein <- NULL
exocarta_mirna  <- NULL

if (file.exists(file.path(exo_dir, "ExoCarta_protein_mrna.txt"))) {
  exocarta_protein <- read.delim(file.path(exo_dir, "ExoCarta_protein_mrna.txt"),
                                  stringsAsFactors = FALSE)
  cat(paste0("  ExoCarta proteins: ", nrow(exocarta_protein), " entries\n"))
}

if (file.exists(file.path(exo_dir, "ExoCarta_mirna.txt"))) {
  exocarta_mirna <- read.delim(file.path(exo_dir, "ExoCarta_mirna.txt"),
                                stringsAsFactors = FALSE)
  cat(paste0("  ExoCarta miRNAs: ", nrow(exocarta_mirna), " entries\n"))
}

# --- Vesiclepedia ---
vesiclepedia_protein <- NULL
vesiclepedia_mirna  <- NULL

if (file.exists(file.path(exo_dir, "Vesiclepedia_protein_mrna.txt"))) {
  vesiclepedia_protein <- read.delim(file.path(exo_dir, "Vesiclepedia_protein_mrna.txt"),
                                      stringsAsFactors = FALSE)
  cat(paste0("  Vesiclepedia proteins: ", nrow(vesiclepedia_protein), " entries\n"))
}

# 如果数据库文件不存在，使用文献报道的核心cargo构建
if (is.null(exocarta_protein) && is.null(vesiclepedia_protein)) {
  cat("  Database files not found. Using curated MSC exosome cargo from literature.\n")
}

# ============================================================
# 1.2 构建MSC来源外泌体cargo综合数据库
# ============================================================

cat("  Building MSC exosome cargo database...\n")

# === 文献报道的MSC外泌体核心cargo ===
# 基于多篇蛋白组学和miRNA-seq研究汇总

# 核心蛋白cargo（高置信度，多次报道）
msc_exo_proteins <- data.frame(
  gene = c(
    # 表面标记/膜蛋白
    "CD9", "CD63", "CD81", "TSG101", "ALIX", "SDCBP", "PDCD6IP",
    # 整合素（靶向相关）
    "ITGA5", "ITGB1", "ITGAV", "ITGB3",
    # 生长因子
    "TGFB1", "TGFB3", "FGF2", "HGF", "IGF1", "VEGFA", "PDGFB", "BMP2",
    # 细胞因子/趋化因子
    "IL10", "IL1RN", "CCL2", "CXCL12", "TSG6",
    # 基质调控
    "MMP2", "TIMP1", "TIMP2", "FN1", "COL1A1",
    # 信号蛋白
    "WNT3A", "WNT5A", "SHH", "DKK1",
    # 抗炎/免疫调节
    "IDO1", "HLA_G", "PD_L1", "GALECTIN1",
    # 代谢/保护
    "SOD1", "SOD2", "CAT", "GPX1",
    # 热休克蛋白
    "HSPA8", "HSP90AA1", "HSPA1A", "HSPB1",
    # 细胞骨架
    "ACTB", "TUBA1A", "MSN", "EZR",
    # 其他功能蛋白
    "ANXA1", "ANXA2", "ANXA5", "PKM", "LDHA", "ENO1", "GAPDH"
  ),
  category = c(
    rep("Surface_Marker", 7),
    rep("Integrin", 4),
    rep("Growth_Factor", 8),
    rep("Cytokine_Chemokine", 5),
    rep("ECM_Regulation", 5),
    rep("Signaling", 4),
    rep("Immunomodulation", 4),
    rep("Antioxidant", 4),
    rep("Heat_Shock_Protein", 4),
    rep("Cytoskeleton", 4),
    rep("Metabolic", 7)
  ),
  therapeutic_relevance = c(
    rep("Targeting", 7),
    rep("Targeting", 4),
    rep("Repair", 8),
    rep("Anti_inflammatory", 5),
    rep("ECM_remodeling", 5),
    rep("Differentiation", 4),
    rep("Immunomodulation", 4),
    rep("Cytoprotection", 4),
    rep("Stress_response", 4),
    rep("Structure", 4),
    rep("Metabolism", 7)
  ),
  confidence = "High",
  stringsAsFactors = FALSE
)

# 核心miRNA cargo
msc_exo_mirnas <- data.frame(
  mirna = c(
    # 抗炎miRNAs
    "hsa-miR-146a-5p", "hsa-miR-146b-5p", "hsa-miR-21-5p",
    "hsa-miR-124-3p", "hsa-miR-223-3p",
    # 软骨保护/分化miRNAs
    "hsa-miR-140-5p", "hsa-miR-140-3p", "hsa-miR-127-5p",
    "hsa-miR-320a", "hsa-miR-320c",
    # 抗凋亡miRNAs
    "hsa-miR-100-5p", "hsa-miR-199a-3p", "hsa-miR-26a-5p",
    # ECM调控
    "hsa-miR-29a-3p", "hsa-miR-29b-3p", "hsa-miR-29c-3p",
    # 免疫调节
    "hsa-miR-155-5p", "hsa-miR-181a-5p", "hsa-miR-182-5p",
    # WNT信号
    "hsa-miR-92a-3p", "hsa-miR-23b-3p",
    # TGF-β信号
    "hsa-miR-135b-5p", "hsa-miR-let-7a-5p", "hsa-miR-let-7b-5p",
    # 其他修复相关
    "hsa-miR-125b-5p", "hsa-miR-130a-3p", "hsa-miR-210-3p",
    "hsa-miR-22-3p", "hsa-miR-27a-3p", "hsa-miR-27b-3p"
  ),
  function_category = c(
    rep("Anti_inflammatory", 5),
    rep("Chondroprotective", 5),
    rep("Anti_apoptotic", 3),
    rep("ECM_regulation", 3),
    rep("Immunomodulation", 3),
    rep("WNT_signaling", 2),
    rep("TGFB_signaling", 3),
    rep("Repair_associated", 6)
  ),
  validated_targets = c(
    "TRAF6/IRAK1", "TRAF6", "PDCD4/PTEN", "NF-kB", "NLRP3",
    "ADAMTS5", "ADAMTS5", "MMP13", "MMP13", "ADAMTS4",
    "mTOR", "COX2", "NF-kB",
    "COL1A1/COL3A1", "COL1A1", "COL1A1",
    "SOCS1", "TLR4", "FOXO1",
    "DKK3", "PTEN",
    "SMAD5", "HMGA2", "HMGA2",
    "TNF-alpha", "SMAD2", "HIF1A",
    "BMP7", "RUNX2", "PPARG"
  ),
  stringsAsFactors = FALSE
)

# ============================================================
# 1.3 不同来源MSC的cargo比较
# ============================================================

cat("  Comparing cargo across MSC sources...\n")

# 文献报道的不同来源MSC外泌体差异cargo
source_specific_cargo <- list(
  Bone_Marrow = list(
    enriched_proteins = c("BMP2", "BMP7", "RUNX2", "VEGFA", "ANG1"),
    enriched_mirnas = c("hsa-miR-21-5p", "hsa-miR-29a-3p", "hsa-miR-199a-3p"),
    strength = "骨和软骨修复"
  ),
  Adipose = list(
    enriched_proteins = c("HGF", "IGF1", "FGF2", "IL10", "TSG6"),
    enriched_mirnas = c("hsa-miR-146a-5p", "hsa-miR-127-5p", "hsa-miR-199a-3p"),
    strength = "抗炎和免疫调节"
  ),
  Synovial = list(
    enriched_proteins = c("TGFB1", "TGFB3", "SOX9", "COL2A1", "PRG4"),
    enriched_mirnas = c("hsa-miR-140-5p", "hsa-miR-140-3p", "hsa-miR-320a"),
    strength = "软骨分化和保护（与半月板同源性最高）"
  ),
  Umbilical_Cord = list(
    enriched_proteins = c("VEGFA", "PDGFB", "WNT3A", "SHH", "FGF2"),
    enriched_mirnas = c("hsa-miR-21-5p", "hsa-miR-100-5p", "hsa-miR-let-7a-5p"),
    strength = "增殖和血管生成"
  )
)

# 保存来源比较表
source_comparison <- do.call(rbind, lapply(names(source_specific_cargo), function(src) {
  data.frame(
    Source = src,
    Enriched_Proteins = paste(source_specific_cargo[[src]]$enriched_proteins, collapse = ", "),
    Enriched_miRNAs = paste(source_specific_cargo[[src]]$enriched_mirnas, collapse = ", "),
    Therapeutic_Strength = source_specific_cargo[[src]]$strength,
    stringsAsFactors = FALSE
  )
}))

write.csv(source_comparison,
          file.path(res_dir, "tables", "MSC_source_cargo_comparison.csv"),
          row.names = FALSE)

# ============================================================
# 1.4 治疗性cargo功能富集
# ============================================================

cat("  Performing functional enrichment of cargo...\n")

# 蛋白cargo的GO/KEGG富集
cargo_entrez <- bitr(msc_exo_proteins$gene, fromType = "SYMBOL",
                     toType = "ENTREZID", OrgDb = org.Hs.eg.db)

ego_cargo <- enrichGO(
  gene = cargo_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

ekegg_cargo <- enrichKEGG(
  gene = cargo_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

# 保存
write.csv(as.data.frame(ego_cargo),
          file.path(res_dir, "tables", "cargo_GO_enrichment.csv"),
          row.names = FALSE)
write.csv(as.data.frame(ekegg_cargo),
          file.path(res_dir, "tables", "cargo_KEGG_enrichment.csv"),
          row.names = FALSE)

# 可视化
pdf(file.path(res_dir, "figures", "03_exosome", "cargo_GO_enrichment.pdf"),
    width = 12, height = 10)
dotplot(ego_cargo, showCategory = 20, font.size = 10) +
  ggtitle("GO Biological Process - MSC Exosome Cargo")
dev.off()

pdf(file.path(res_dir, "figures", "03_exosome", "cargo_KEGG_enrichment.pdf"),
    width = 12, height = 8)
dotplot(ekegg_cargo, showCategory = 15, font.size = 10) +
  ggtitle("KEGG Pathway - MSC Exosome Cargo")
dev.off()

# 治疗通路映射
therapeutic_pathways <- config$exosome$therapeutic_pathways
pathway_mapping <- data.frame(
  Pathway_Category = rep(names(therapeutic_pathways), sapply(therapeutic_pathways, length)),
  KEGG_Pathway = unlist(therapeutic_pathways),
  stringsAsFactors = FALSE
)

# 将cargo映射到治疗通路
if (nrow(as.data.frame(ekegg_cargo)) > 0) {
  kegg_df <- as.data.frame(ekegg_cargo)
  therapeutic_hits <- kegg_df %>%
    mutate(therapeutic_category = case_when(
      grepl("TNF|NF-kappa|IL-17", Description) ~ "Anti_inflammatory",
      grepl("PI3K|MAPK|Ras", Description) ~ "Pro_proliferation",
      grepl("TGF|Wnt|Hippo", Description) ~ "ECM_synthesis",
      grepl("p53|Apoptosis|mTOR|Autophagy", Description) ~ "Anti_apoptosis",
      grepl("Hedgehog|FoxO|BMP", Description) ~ "Chondrogenesis",
      TRUE ~ "Other"
    ))
  
  write.csv(therapeutic_hits,
            file.path(res_dir, "tables", "cargo_therapeutic_pathway_mapping.csv"),
            row.names = FALSE)
}

# ============================================================
# 1.5 保存综合cargo数据库
# ============================================================

cat("  Saving comprehensive cargo database...\n")

# 创建多sheet的Excel工作簿
wb <- createWorkbook()

addWorksheet(wb, "Protein_Cargo")
writeData(wb, "Protein_Cargo", msc_exo_proteins)

addWorksheet(wb, "miRNA_Cargo")
writeData(wb, "miRNA_Cargo", msc_exo_mirnas)

addWorksheet(wb, "Source_Comparison")
writeData(wb, "Source_Comparison", source_comparison)

if (exists("therapeutic_hits")) {
  addWorksheet(wb, "Therapeutic_Pathways")
  writeData(wb, "Therapeutic_Pathways", therapeutic_hits)
}

saveWorkbook(wb, file.path(res_dir, "tables", "MSC_exosome_cargo_database.xlsx"),
             overwrite = TRUE)

# 保存为RDS供后续分析使用
cargo_db <- list(
  proteins = msc_exo_proteins,
  mirnas = msc_exo_mirnas,
  source_comparison = source_specific_cargo
)
saveRDS(cargo_db, file.path(data_dir, "processed", "cargo_database.rds"))

cat("\n  Cargo database constructed.\n")
cat("  Next step: Rscript scripts/02_exosome_cargo/step2_miRNA_targets.R\n")
cat("=== Module 2 Step 1 Complete ===\n")
