options(conflicts.policy = list(warn = FALSE))
suppressPackageStartupMessages({
  library(Seurat); library(CellChat); library(yaml); library(ggplot2)
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(ComplexHeatmap); library(circlize); library(igraph)
})

config <- yaml::read_yaml("config/params.yaml")
data_dir <- config$data_dir; res_dir <- config$results_dir
merged <- readRDS(file.path(data_dir, "processed/meniscus_annotated.rds"))
cargo_db <- readRDS(file.path(data_dir, "processed/cargo_database.rds"))
mirna_analysis <- readRDS(file.path(data_dir, "processed/mirna_target_analysis.rds"))

cat("=== Rerunning Module 3 with corrected conditions ===\n")
cat("Conditions:", paste(names(table(merged$condition)), collapse=", "), "\n\n")

cat("Part 1: Receptor fingerprinting...\n")
CellChatDB <- CellChatDB.human
all_receptors <- unique(c(
  unlist(strsplit(CellChatDB$interaction$receptor, "_")),
  "TGFBR1","TGFBR2","BMPR1A","FGFR1","FGFR2","FGFR3",
  "IL1R1","IL6R","TNFRSF1A","ITGA5","ITGB1","CD44",
  "PDGFRA","PDGFRB","KDR","MET","EGFR","NOTCH1","FZD1","FZD2","PTCH1"))
expressed_receptors <- all_receptors[all_receptors %in% rownames(merged)]
cell_types <- unique(merged$cell_type)
receptor_expr <- as.matrix(AverageExpression(merged, features=expressed_receptors, group.by="cell_type")$RNA)

cat("  Differential receptors per cell type...\n")
deg_receptors <- list()
for (ct in cell_types) {
  ct_cells <- subset(merged, cell_type == ct)
  if (min(table(ct_cells$condition)) < 10) next
  Idents(ct_cells) <- "condition"
  tryCatch({
    degs <- FindMarkers(ct_cells, ident.1="Degenerated", ident.2="Normal",
                        features=expressed_receptors, min.pct=0.05, logfc.threshold=0.1, verbose=FALSE)
    if (nrow(degs) > 0) { degs$receptor <- rownames(degs); degs$cell_type <- ct; deg_receptors[[ct]] <- degs }
  }, error = function(e) NULL)
}
if (length(deg_receptors) > 0) {
  deg_receptors_df <- bind_rows(deg_receptors)
  write.csv(deg_receptors_df, file.path(res_dir, "tables/differentially_expressed_receptors.csv"), row.names=FALSE)
  cat(paste0("  DE receptors found in ", length(deg_receptors), " cell types\n"))
}

receptor_zscore <- t(scale(t(receptor_expr))); receptor_zscore[is.na(receptor_zscore)] <- 0
pdf(file.path(res_dir, "figures/04_docking/receptor_fingerprint_heatmap.pdf"), width=20, height=18)
Heatmap(receptor_zscore, name="Z-score",
  col=colorRamp2(c(-2,0,2), c("#4DBBD5","white","#E64B35")),
  show_row_names=TRUE, show_column_names=TRUE, column_names_rot=45,
  row_names_gp=gpar(fontsize=4), column_title="Cell Type-Specific Receptor Fingerprint", row_title="Receptors")
dev.off()

cat("\nPart 2: Cargo-receptor matching...\n")
cargo_proteins <- cargo_db$proteins$gene
interaction_db <- CellChatDB$interaction
cargo_ligand_matches <- interaction_db[sapply(strsplit(interaction_db$ligand,"_"), function(x) any(x %in% cargo_proteins)),]

protein_scores <- matrix(0, nrow=length(cargo_proteins), ncol=length(cell_types),
                         dimnames=list(cargo_proteins, cell_types))
for (i in seq_len(nrow(cargo_ligand_matches))) {
  ligand_genes <- unlist(strsplit(cargo_ligand_matches$ligand[i], "_"))
  receptor_genes <- unlist(strsplit(cargo_ligand_matches$receptor[i], "_"))
  matched_cargo <- intersect(ligand_genes, cargo_proteins)
  matched_receptors <- intersect(receptor_genes, rownames(receptor_expr))
  if (length(matched_cargo) > 0 && length(matched_receptors) > 0) {
    for (ct in cell_types) {
      if (ct %in% colnames(receptor_expr)) {
        rs <- mean(receptor_expr[matched_receptors, ct, drop=FALSE], na.rm=TRUE)
        for (cargo in matched_cargo) protein_scores[cargo, ct] <- protein_scores[cargo, ct] + rs
      }
    }
  }
}

mirna_targets <- mirna_analysis$high_conf_targets
target_genes <- unique(mirna_targets$Target)
target_genes <- target_genes[target_genes %in% rownames(merged)]
target_expr <- as.matrix(AverageExpression(merged, features=target_genes, group.by="cell_type")$RNA)
mirna_ids <- unique(mirna_targets$miRNA)
mirna_scores <- matrix(0, nrow=length(mirna_ids), ncol=length(cell_types), dimnames=list(mirna_ids, cell_types))
for (mir in mirna_ids) {
  tgts <- mirna_targets$Target[mirna_targets$miRNA == mir]
  tgts <- tgts[tgts %in% rownames(target_expr)]
  if (length(tgts) > 0) for (ct in cell_types) if (ct %in% colnames(target_expr))
    mirna_scores[mir, ct] <- mean(target_expr[tgts, ct, drop=FALSE], na.rm=TRUE)
}

normalize_mat <- function(m) { r <- range(m, na.rm=TRUE); if(r[2]==r[1]) m*0 else (m-r[1])/(r[2]-r[1]) }
protein_scores_norm <- normalize_mat(protein_scores)
mirna_scores_norm <- normalize_mat(mirna_scores)
write.csv(protein_scores_norm, file.path(res_dir, "tables/protein_cargo_celltype_scores.csv"))
write.csv(mirna_scores_norm, file.path(res_dir, "tables/miRNA_cargo_celltype_scores.csv"))

cat("\nPart 3: CellChat per condition...\n")
run_cellchat <- function(seurat_obj, label) {
  data_input <- GetAssayData(seurat_obj, layer="data")
  meta <- seurat_obj@meta.data; meta$labels <- meta$cell_type
  cc <- createCellChat(object=data_input, meta=meta, group.by="labels")
  cc@DB <- subsetDB(CellChatDB.human, search=c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"))
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type="triMean")
  cc <- filterCommunication(cc, min.cells=10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cat(paste0("  ", label, ": ", length(cc@netP$pathways), " pathways\n"))
  return(cc)
}

cellchat_all <- run_cellchat(merged, "All")
pdf(file.path(res_dir, "figures/04_docking/cellchat_network.pdf"), width=12, height=10)
netVisual_circle(cellchat_all@net$count, vertex.weight=table(cellchat_all@idents),
  weight.scale=TRUE, label.edge=FALSE, title.name="Communication Network")
dev.off()

normal_obj <- subset(merged, condition=="Normal")
degen_obj <- subset(merged, condition=="Degenerated")
cellchat_normal <- tryCatch(run_cellchat(normal_obj, "Normal"), error=function(e){cat("  Normal CellChat failed\n"); NULL})
cellchat_degen <- tryCatch(run_cellchat(degen_obj, "Degenerated"), error=function(e){cat("  Degen CellChat failed\n"); NULL})

lr_pairs <- cellchat_all@LR$LRsig
exo_repair <- data.frame(Pathway=character(), Cargo=character(), Receptor=character(), stringsAsFactors=FALSE)
for (i in seq_len(nrow(lr_pairs))) {
  ligands <- unlist(strsplit(lr_pairs$ligand[i], "_"))
  cargo_match <- intersect(ligands, cargo_proteins)
  if (length(cargo_match) > 0)
    exo_repair <- rbind(exo_repair, data.frame(Pathway=lr_pairs$pathway_name[i],
      Cargo=paste(cargo_match,collapse=","), Receptor=lr_pairs$receptor[i]))
}
exo_repair <- distinct(exo_repair)
write.csv(exo_repair, file.path(res_dir, "tables/exosome_repair_map.csv"), row.names=FALSE)

cat("\nPart 4: Treatability scores...\n")
treat <- data.frame(cell_type=cell_types,
  protein_match=colMeans(protein_scores_norm, na.rm=TRUE),
  mirna_target=colMeans(mirna_scores_norm, na.rm=TRUE))
for (col in c("protein_match","mirna_target")) {
  rng <- range(treat[[col]], na.rm=TRUE)
  if (rng[2]>rng[1]) treat[[col]] <- (treat[[col]]-rng[1])/(rng[2]-rng[1])
}
treat$composite <- rowMeans(treat[,c("protein_match","mirna_target")], na.rm=TRUE)
treat <- treat %>% arrange(desc(composite))
write.csv(treat, file.path(res_dir, "tables/treatability_scores.csv"), row.names=FALSE)

docking_results <- list(receptor_expr=receptor_expr, protein_scores=protein_scores_norm,
  mirna_scores=mirna_scores_norm, cellchat_all=cellchat_all,
  cellchat_normal=cellchat_normal, cellchat_degen=cellchat_degen,
  exo_repair_map=exo_repair, treatability_scores=treat)
saveRDS(docking_results, file.path(data_dir, "processed/docking_results.rds"))

summary_df <- data.frame(
  Metric=c("Total cells","Normal cells","Degenerated cells","Clusters","Active pathways",
           "Exosome repair pairs","Top treatable cluster"),
  Value=c(ncol(merged), sum(merged$condition=="Normal"), sum(merged$condition=="Degenerated"),
          length(cell_types), length(cellchat_all@netP$pathways), nrow(exo_repair), treat$cell_type[1]))
write.csv(summary_df, file.path(res_dir, "tables/analysis_summary.csv"), row.names=FALSE)
print(summary_df)

cat("\n=== Module 3 Updated Complete ===\n")
