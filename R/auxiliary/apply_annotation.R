# 正式注释脚本 — 下次AutoDL开机后执行
# Rscript /root/autodl-tmp/meniscus_exosome_project/scripts/apply_annotation.R

library(Seurat); library(ggplot2); library(yaml); library(dplyr)
config <- yaml::read_yaml("config/params.yaml")
data_dir <- config$data_dir; res_dir <- config$results_dir

obj <- readRDS(file.path(data_dir, "processed/meniscus_annotated.rds"))

# 正式注释
annotation_map <- c(
  "0"  = "FC-Regulatory",
  "1"  = "FC-ECM",
  "2"  = "Secretory_Chondrocyte",
  "3"  = "FC-Metabolic",
  "4"  = "Meniscus_Progenitor",
  "5"  = "FC-Degenerated",
  "6"  = "Proliferating_G2M",
  "7"  = "FC-Inflammatory",
  "8"  = "Proliferating_M",
  "9"  = "Proliferating_S",
  "10" = "Male_Inflammatory",
  "11" = "Smooth_Muscle",
  "12" = "Hypertrophic_Chondrocyte",
  "13" = "OuterZone_FC",
  "14" = "Macrophage",
  "15" = "Endothelial",
  "16" = "Stem_Progenitor"
)

obj@meta.data$cell_type <- annotation_map[as.character(obj$seurat_clusters)]

# 大类分组
obj@meta.data$cell_category <- case_when(
  grepl("FC-|Chondrocyte|OuterZone", obj$cell_type) ~ "Fibrochondrocyte",
  grepl("Proliferating", obj$cell_type)             ~ "Proliferating",
  grepl("Progenitor|Stem", obj$cell_type)           ~ "Progenitor/Stem",
  grepl("Macrophage|Inflammatory", obj$cell_type)   ~ "Immune/Inflammatory",
  TRUE                                              ~ "Stromal/Vascular"
)

cat("=== Cell Type Composition ===\n")
print(sort(table(obj$cell_type), decreasing=TRUE))

cat("\n=== Category Composition ===\n")
print(sort(table(obj$cell_category), decreasing=TRUE))

# 保存
saveRDS(obj, file.path(data_dir, "processed/meniscus_annotated.rds"))
saveRDS(obj, file.path(data_dir, "processed/meniscus_integrated.rds"))

# --- 论文级 UMAP ---
colors17 <- c(
  "FC-Regulatory"="#E64B35", "FC-ECM"="#F39B7F", "Secretory_Chondrocyte"="#D95F02",
  "FC-Metabolic"="#E7298A", "Meniscus_Progenitor"="#1B9E77", "FC-Degenerated"="#7570B3",
  "Proliferating_G2M"="#66A61E", "FC-Inflammatory"="#E6AB02", "Proliferating_M"="#A6D854",
  "Proliferating_S"="#B3DE69", "Male_Inflammatory"="#BC80BD", "Smooth_Muscle"="#8DD3C7",
  "Hypertrophic_Chondrocyte"="#FB8072", "OuterZone_FC"="#FDB462",
  "Macrophage"="#4DBBD5", "Endothelial"="#00A087", "Stem_Progenitor"="#3C5488"
)

pdf(file.path(res_dir, "figures/01_integration/UMAP_annotated.pdf"), width=14, height=10)
print(DimPlot(obj, group.by="cell_type", label=TRUE, repel=TRUE,
              cols=colors17, label.size=3.5, pt.size=0.3) +
      ggtitle("Human Meniscus Single-Cell Atlas") +
      theme(legend.text=element_text(size=9)))
dev.off()

pdf(file.path(res_dir, "figures/01_integration/UMAP_category.pdf"), width=12, height=8)
print(DimPlot(obj, group.by="cell_category", label=TRUE, repel=TRUE, pt.size=0.3) +
      ggtitle("Cell Category Overview"))
dev.off()

pdf(file.path(res_dir, "figures/01_integration/UMAP_condition_annotated.pdf"), width=24, height=10)
print(DimPlot(obj, group.by="cell_type", split.by="condition", label=TRUE, repel=TRUE,
              cols=colors17, label.size=3, pt.size=0.3) +
      ggtitle("Normal vs Degenerated"))
dev.off()

# 细胞比例变化
prop <- obj@meta.data %>%
  group_by(condition, cell_type) %>% summarise(n=n(), .groups="drop") %>%
  group_by(condition) %>% mutate(prop=n/sum(n))

pdf(file.path(res_dir, "figures/01_integration/proportion_barplot.pdf"), width=14, height=6)
print(ggplot(prop, aes(x=cell_type, y=prop, fill=condition)) +
      geom_bar(stat="identity", position="dodge") +
      theme_classic() + coord_flip() +
      labs(title="Cell Type Proportion: Normal vs Degenerated", x="", y="Proportion") +
      scale_fill_manual(values=c("Normal"="#4DBBD5","Degenerated"="#E64B35")))
dev.off()

cat("\n=== Annotation Complete ===\n")
cat("Figures saved to results/figures/01_integration/\n")
