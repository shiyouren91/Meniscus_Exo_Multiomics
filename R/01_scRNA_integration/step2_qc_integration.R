# ============================================================
# Step 2: 单细胞数据质控与整合
# 
# 输入: 原始count matrix
# 输出: 整合后的Seurat对象 + QC报告
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(DoubletFinder)
  library(yaml)
  library(tidyverse)
  library(patchwork)
  library(future)
})

# 配置并行
plan("multicore", workers = 4)
options(future.globals.maxSize = 50 * 1024^3)  # 50GB

# 读取配置
config <- yaml::read_yaml("config/params.yaml")
data_dir  <- path.expand(config$data_dir)
res_dir   <- path.expand(config$results_dir)
qc_params <- config$scRNA$qc
int_params <- config$scRNA$integration

dir.create(file.path(res_dir, "figures", "01_integration"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(res_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

cat("=== Step 2: QC & Integration ===\n")

# ============================================================
# 2.1 加载各数据集
# ============================================================

load_dataset <- function(dataset_id, data_path) {
  cat(paste0("  Loading ", dataset_id, "...\n"))
  
  # 尝试不同的数据格式
  if (dir.exists(file.path(data_path, "filtered_feature_bc_matrix"))) {
    # 10X standard output
    mat <- Read10X(file.path(data_path, "filtered_feature_bc_matrix"))
  } else if (dir.exists(file.path(data_path, "outs"))) {
    mat <- Read10X(file.path(data_path, "outs", "filtered_feature_bc_matrix"))
  } else {
    # 尝试读取 .h5 文件
    h5_files <- list.files(data_path, pattern = "\\.h5$", full.names = TRUE)
    if (length(h5_files) > 0) {
      mat <- Read10X_h5(h5_files[1])
    } else {
      # 尝试读取 matrix.mtx + genes.tsv + barcodes.tsv
      mtx_files <- list.files(data_path, pattern = "matrix", full.names = TRUE, recursive = TRUE)
      if (length(mtx_files) > 0) {
        mtx_dir <- dirname(mtx_files[1])
        mat <- Read10X(mtx_dir)
      } else {
        # 尝试读取 RDS 文件
        rds_files <- list.files(data_path, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
        if (length(rds_files) > 0) {
          return(readRDS(rds_files[1]))
        }
        stop(paste("Cannot find data files in", data_path))
      }
    }
  }
  
  obj <- CreateSeuratObject(
    counts = mat,
    project = dataset_id,
    min.cells = qc_params$min_cells,
    min.features = qc_params$min_genes
  )
  obj$dataset <- dataset_id
  return(obj)
}

# 加载所有可用数据集
seurat_list <- list()
dataset_dirs <- list.dirs(file.path(data_dir, "raw", "scRNA"), recursive = FALSE)

for (d in dataset_dirs) {
  dataset_id <- basename(d)
  # 跳过空目录
  if (length(list.files(d, recursive = TRUE)) == 0) {
    cat(paste0("  Skipping empty directory: ", dataset_id, "\n"))
    next
  }
  
  # 如果有多个样本子目录（如GSE133449解压后多个样本）
  sample_dirs <- list.dirs(d, recursive = FALSE)
  if (length(sample_dirs) > 0 && all(sapply(sample_dirs, function(x) length(list.files(x)) > 0))) {
    for (sd in sample_dirs) {
      sample_id <- paste0(dataset_id, "_", basename(sd))
      tryCatch({
        seurat_list[[sample_id]] <- load_dataset(sample_id, sd)
      }, error = function(e) {
        cat(paste0("  Warning: Failed to load ", sample_id, ": ", e$message, "\n"))
      })
    }
  } else {
    tryCatch({
      seurat_list[[dataset_id]] <- load_dataset(dataset_id, d)
    }, error = function(e) {
      cat(paste0("  Warning: Failed to load ", dataset_id, ": ", e$message, "\n"))
    })
  }
}

cat(paste0("  Loaded ", length(seurat_list), " datasets\n"))

# ============================================================
# 2.2 质控
# ============================================================

cat("\n=== Quality Control ===\n")

qc_one_sample <- function(obj) {
  obj[["percent.mt"]]  <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  
  # 记录质控前细胞数
  n_before <- ncol(obj)
  
  obj <- subset(obj,
    nFeature_RNA > qc_params$min_genes &
    nFeature_RNA < qc_params$max_genes &
    percent.mt < qc_params$max_percent_mt
  )
  
  n_after <- ncol(obj)
  cat(paste0("    ", obj$dataset[1], ": ", n_before, " -> ", n_after, 
             " cells (removed ", n_before - n_after, ")\n"))
  return(obj)
}

seurat_list <- lapply(seurat_list, qc_one_sample)

# 质控可视化
qc_plots <- lapply(seurat_list, function(obj) {
  VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
          ncol = 3, pt.size = 0) + 
    ggtitle(obj$dataset[1])
})

pdf(file.path(res_dir, "figures", "01_integration", "QC_violin_plots.pdf"),
    width = 14, height = 4 * length(qc_plots))
wrap_plots(qc_plots, ncol = 1)
dev.off()

# ============================================================
# 2.3 DoubletFinder 去双联体
# ============================================================

cat("\n=== Doublet Detection ===\n")

remove_doublets <- function(obj) {
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  
  # DoubletFinder
  sweep_res <- paramSweep(obj, PCs = 1:20, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  
  optimal_pK <- as.numeric(as.character(
    bcmvn$pK[which.max(bcmvn$BCmetric)]
  ))
  
  nExp_poi <- round(qc_params$doublet_rate * ncol(obj))
  
  obj <- doubletFinder(obj, PCs = 1:20, pN = 0.25, pK = optimal_pK,
                       nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  # 获取 DoubletFinder 列名
  df_col <- grep("^DF.classifications", colnames(obj@meta.data), value = TRUE)[1]
  
  n_doublets <- sum(obj@meta.data[[df_col]] == "Doublet")
  cat(paste0("    ", obj$dataset[1], ": removed ", n_doublets, " doublets\n"))
  
  obj <- subset(obj, cells = colnames(obj)[obj@meta.data[[df_col]] == "Singlet"])
  
  # 清理中间列
  obj@meta.data <- obj@meta.data[, !grepl("^pANN|^DF\\.", colnames(obj@meta.data))]
  
  return(obj)
}

seurat_list <- lapply(seurat_list, function(obj) {
  tryCatch(remove_doublets(obj), error = function(e) {
    cat(paste0("    Warning: DoubletFinder failed for ", obj$dataset[1], ", skipping\n"))
    return(obj)
  })
})

# ============================================================
# 2.4 合并与整合
# ============================================================

cat("\n=== Merging & Integration ===\n")

# 合并所有对象
if (length(seurat_list) > 1) {
  merged <- merge(seurat_list[[1]], y = seurat_list[-1],
                  add.cell.ids = names(seurat_list))
} else {
  merged <- seurat_list[[1]]
}

cat(paste0("  Total cells after merge: ", ncol(merged), "\n"))
cat(paste0("  Total genes: ", nrow(merged), "\n"))

# 标准预处理
merged <- NormalizeData(merged, verbose = FALSE)
merged <- FindVariableFeatures(merged, selection.method = "vst",
                               nfeatures = 3000, verbose = FALSE)
merged <- ScaleData(merged, vars.to.regress = c("percent.mt", "nCount_RNA"),
                    verbose = FALSE)
merged <- RunPCA(merged, npcs = 50, verbose = FALSE)

# ElbowPlot
pdf(file.path(res_dir, "figures", "01_integration", "ElbowPlot.pdf"), width = 8, height = 5)
ElbowPlot(merged, ndims = 50)
dev.off()

# Harmony 整合
cat("  Running Harmony integration...\n")
merged <- RunHarmony(merged,
                     group.by.vars = "dataset",
                     theta = int_params$harmony_theta,
                     lambda = int_params$harmony_lambda,
                     max.iter.harmony = 30,
                     verbose = FALSE)

# 降维与聚类
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:int_params$n_pcs, verbose = FALSE)
merged <- RunTSNE(merged, reduction = "harmony", dims = 1:int_params$n_pcs, verbose = FALSE)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:int_params$n_pcs, verbose = FALSE)

# 多分辨率聚类
for (res in config$scRNA$clustering$resolution) {
  merged <- FindClusters(merged, resolution = res, algorithm = 4, verbose = FALSE)  # 4 = Leiden
}

# 设定默认分辨率
default_res <- config$scRNA$clustering$selected_resolution
Idents(merged) <- paste0("RNA_snn_res.", default_res)

# ============================================================
# 2.5 整合效果可视化
# ============================================================

cat("\n=== Generating integration QC plots ===\n")

# UMAP by cluster and dataset
p1 <- DimPlot(merged, reduction = "umap", label = TRUE, repel = TRUE,
              label.size = 4) + NoLegend() + ggtitle("Clusters")
p2 <- DimPlot(merged, reduction = "umap", group.by = "dataset") +
  ggtitle("By Dataset")

pdf(file.path(res_dir, "figures", "01_integration", "UMAP_integration.pdf"),
    width = 16, height = 7)
p1 + p2
dev.off()

# 各数据集在每个cluster中的比例
prop_table <- table(Idents(merged), merged$dataset) %>%
  as.data.frame() %>%
  rename(Cluster = Var1, Dataset = Var2, Count = Freq) %>%
  group_by(Cluster) %>%
  mutate(Proportion = Count / sum(Count))

p3 <- ggplot(prop_table, aes(x = Cluster, y = Proportion, fill = Dataset)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic() +
  labs(title = "Dataset Proportion per Cluster", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(res_dir, "figures", "01_integration", "dataset_proportion.pdf"),
    width = 12, height = 5)
print(p3)
dev.off()

# ============================================================
# 2.6 保存结果
# ============================================================

cat("\n=== Saving results ===\n")

# 保存整合后Seurat对象
saveRDS(merged, file.path(data_dir, "processed", "meniscus_integrated.rds"))

# 保存metadata
write.csv(merged@meta.data,
          file.path(res_dir, "tables", "cell_metadata_integrated.csv"))

# 保存细胞统计
cell_stats <- data.frame(
  Dataset = names(table(merged$dataset)),
  Cells = as.numeric(table(merged$dataset))
)
write.csv(cell_stats, file.path(res_dir, "tables", "cell_counts_by_dataset.csv"),
          row.names = FALSE)

cat(paste0("\n  Integrated object saved: ", ncol(merged), " cells, ", nrow(merged), " genes\n"))
cat("  Next step: Rscript scripts/01_scRNA_integration/step3_annotation.R\n")
cat("=== Step 2 Complete ===\n")
