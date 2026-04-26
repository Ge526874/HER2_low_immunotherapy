# Peer-review R script: immune-cell single-cell RNA-seq analysis

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(readxl)
  library(openxlsx)
  library(DESeq2)
  library(EnhancedVolcano)
  library(msigdbr)
  library(clusterProfiler)
})

set.seed(123)

input_seurat_file <- "input/seurat_object.rds"
input_clinical_file <- "input/clinical_metadata.xlsx"
output_dir <- "peer_review_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

seurat_object <- readRDS(input_seurat_file)
clinical_metadata <- read_excel(input_clinical_file) |> as.data.frame()

required_clinical_columns <- c(
  "sample_id", "patient_id", "HER2_group", "timepoint",
  "treatment_group", "response_group", "ER_status"
)

missing_clinical_columns <- setdiff(required_clinical_columns, colnames(clinical_metadata))
if (length(missing_clinical_columns) > 0) {
  stop("Missing required clinical columns: ", paste(missing_clinical_columns, collapse = ", "))
}

rownames(clinical_metadata) <- clinical_metadata$sample_id

seurat_object$patient_id <- clinical_metadata[seurat_object$sample_id, "patient_id"]
seurat_object$HER2_group <- clinical_metadata[seurat_object$sample_id, "HER2_group"]
seurat_object$timepoint <- clinical_metadata[seurat_object$sample_id, "timepoint"]
seurat_object$treatment_group <- clinical_metadata[seurat_object$sample_id, "treatment_group"]
seurat_object$response_group <- clinical_metadata[seurat_object$sample_id, "response_group"]
seurat_object$ER_status <- clinical_metadata[seurat_object$sample_id, "ER_status"]

immune_major_types <- c("T cell", "NK cell", "B cell", "Plasma cell", "Myeloid", "DC", "Mast")
immune_object <- subset(seurat_object, subset = major_celltype %in% immune_major_types)

DefaultAssay(immune_object) <- "RNA"

if (!"umap" %in% names(immune_object@reductions)) {
  immune_object <- NormalizeData(immune_object)
  immune_object <- FindVariableFeatures(immune_object)
  immune_object <- ScaleData(immune_object)
  immune_object <- RunPCA(immune_object)
  immune_object <- FindNeighbors(immune_object, dims = 1:30)
  immune_object <- FindClusters(immune_object, resolution = 1)
  immune_object <- RunUMAP(immune_object, dims = 1:30)
}

p_umap_major <- DimPlot(
  immune_object,
  reduction = "umap",
  group.by = "major_celltype",
  label = TRUE,
  raster = FALSE
)

ggsave(
  filename = file.path(output_dir, "immune_major_celltype_umap.pdf"),
  plot = p_umap_major,
  width = 8,
  height = 6
)

immune_marker_list <- list(
  "CD4_T" = c("CD3D", "CD3E", "CD4", "IL7R", "TCF7"),
  "CD8_T" = c("CD3D", "CD3E", "CD8A", "CD8B", "GZMB", "NKG7"),
  "Treg" = c("FOXP3", "IL2RA", "CTLA4"),
  "NK" = c("NKG7", "GNLY", "KLRD1", "FCGR3A"),
  "B_cell" = c("MS4A1", "CD79A", "CD79B"),
  "Plasma_cell" = c("MZB1", "JCHAIN", "XBP1"),
  "Myeloid" = c("LYZ", "S100A8", "S100A9", "FCGR3A"),
  "Macrophage" = c("CD68", "CD163", "MRC1"),
  "DC" = c("FCER1A", "CLEC10A", "LILRA4"),
  "Mast" = c("TPSAB1", "TPSB2", "CPA3")
)

available_marker_list <- lapply(
  immune_marker_list,
  function(x) intersect(x, rownames(immune_object))
)
available_marker_list <- available_marker_list[lengths(available_marker_list) > 0]

p_marker_dotplot <- DotPlot(
  immune_object,
  features = available_marker_list,
  group.by = "major_celltype"
) +
  RotatedAxis()

ggsave(
  filename = file.path(output_dir, "immune_marker_dotplot.pdf"),
  plot = p_marker_dotplot,
  width = 12,
  height = 6
)

cell_fraction <- immune_object@meta.data |>
  filter(!is.na(HER2_group), !is.na(timepoint), !is.na(patient_id)) |>
  count(patient_id, sample_id, HER2_group, timepoint, treatment_group, response_group, major_celltype, name = "cell_count") |>
  group_by(patient_id, sample_id, HER2_group, timepoint, treatment_group, response_group) |>
  mutate(fraction = cell_count / sum(cell_count)) |>
  ungroup()

write.xlsx(cell_fraction, file.path(output_dir, "immune_cell_fraction_by_sample.xlsx"), rowNames = FALSE)

p_fraction <- ggplot(
  cell_fraction |> filter(timepoint == "pre"),
  aes(x = major_celltype, y = fraction, fill = HER2_group)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(output_dir, "pretreatment_immune_cell_fraction.pdf"),
  plot = p_fraction,
  width = 10,
  height = 5
)

paired_fraction_change <- cell_fraction |>
  select(patient_id, HER2_group, treatment_group, response_group, major_celltype, timepoint, fraction) |>
  distinct() |>
  pivot_wider(names_from = timepoint, values_from = fraction) |>
  filter(!is.na(pre), !is.na(post)) |>
  mutate(fraction_change = post - pre)

write.xlsx(paired_fraction_change, file.path(output_dir, "paired_immune_fraction_change.xlsx"), rowNames = FALSE)

p_fraction_change <- ggplot(
  paired_fraction_change,
  aes(x = major_celltype, y = fraction_change, fill = HER2_group)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(output_dir, "paired_immune_fraction_change.pdf"),
  plot = p_fraction_change,
  width = 10,
  height = 5
)

immune_gene_sets <- list(
  "Antigen_presentation" = c("HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "B2M", "TAP1", "TAP2"),
  "Immune_checkpoint" = c("CD274", "PDCD1", "PDCD1LG2", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "IDO1"),
  "Cytotoxicity" = c("GZMB", "GZMA", "PRF1", "NKG7", "GNLY", "IFNG")
)

run_seurat_de <- function(object, celltypes, group_label, ident_1, ident_2, output_prefix) {
  subset_object <- subset(object, subset = major_celltype %in% celltypes)
  subset_object$comparison_group <- paste(group_label, subset_object$HER2_group, subset_object$timepoint, sep = "_")
  Idents(subset_object) <- "comparison_group"

  de_result <- FindMarkers(
    subset_object,
    ident.1 = ident_1,
    ident.2 = ident_2,
    logfc.threshold = 0,
    min.pct = 0.05
  )

  de_result$gene <- rownames(de_result)
  write.xlsx(de_result, file.path(output_dir, paste0(output_prefix, "_seurat_de.xlsx")), rowNames = FALSE)

  de_result$p_val_adj[de_result$p_val_adj == 0] <- min(de_result$p_val_adj[de_result$p_val_adj > 0], na.rm = TRUE)

  p_volcano <- EnhancedVolcano(
    de_result,
    lab = de_result$gene,
    x = "avg_log2FC",
    y = "p_val_adj",
    pCutoff = 0.05,
    FCcutoff = 0.25
  )

  ggsave(
    filename = file.path(output_dir, paste0(output_prefix, "_seurat_de_volcano.pdf")),
    plot = p_volcano,
    width = 8,
    height = 8
  )

  return(de_result)
}

run_pseudobulk_de <- function(object, celltypes, output_prefix) {
  subset_object <- subset(object, subset = major_celltype %in% celltypes)

  pseudobulk_object <- AggregateExpression(
    subset_object,
    assays = "RNA",
    return.seurat = TRUE,
    group.by = c("HER2_group", "timepoint", "patient_id")
  )

  count_matrix <- GetAssayData(pseudobulk_object, assay = "RNA", slot = "counts")
  sample_metadata <- pseudobulk_object@meta.data

  sample_metadata$HER2_timepoint <- paste(sample_metadata$HER2_group, sample_metadata$timepoint, sep = "_")
  sample_metadata$HER2_timepoint <- factor(sample_metadata$HER2_timepoint)
  sample_metadata$patient_id <- factor(sample_metadata$patient_id)

  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(count_matrix)),
    colData = sample_metadata,
    design = ~ patient_id + HER2_timepoint
  )

  dds <- DESeq(dds)
  comparison_levels <- levels(sample_metadata$HER2_timepoint)

  if (!all(c("1+_pre", "0_pre") %in% comparison_levels)) {
    stop("Expected comparison groups '1+_pre' and '0_pre' were not found.")
  }

  de_result <- results(dds, contrast = c("HER2_timepoint", "1+_pre", "0_pre"))
  de_result <- as.data.frame(de_result)
  de_result$gene <- rownames(de_result)

  write.xlsx(de_result, file.path(output_dir, paste0(output_prefix, "_pseudobulk_de.xlsx")), rowNames = FALSE)

  de_result$padj[de_result$padj == 0] <- min(de_result$padj[de_result$padj > 0], na.rm = TRUE)

  p_volcano <- EnhancedVolcano(
    de_result,
    lab = de_result$gene,
    x = "log2FoldChange",
    y = "padj",
    pCutoff = 0.05,
    FCcutoff = 0.25
  )

  ggsave(
    filename = file.path(output_dir, paste0(output_prefix, "_pseudobulk_de_volcano.pdf")),
    plot = p_volcano,
    width = 8,
    height = 8
  )

  return(de_result)
}

celltype_groups <- list(
  "cd4_t" = c("CD4+T", "CD4+TCF1+Tn", "CD4+PD1+Tex", "CD4+TCF1+Tem", "CD4+Treg"),
  "cd8_t" = c("CD8+T", "CD8+TCF1+Tem", "CD8+GZMB+PD1+Teff", "CD8+Trm", "CD8+Tprf"),
  "nk" = c("CD56+NK"),
  "b_cell" = c("CD79A+pB", "CD20+Bn", "Bprf"),
  "apc_myeloid" = c("PD.L1+IDO+APCs", "M2.Macro", "M1.Macro", "DCs", "Macro.prf", "M2.TAM")
)

if ("immune_subtype" %in% colnames(immune_object@meta.data)) {
  immune_object$major_celltype <- immune_object$immune_subtype
}

de_results <- list()
pseudobulk_results <- list()

for (analysis_name in names(celltype_groups)) {
  current_celltypes <- celltype_groups[[analysis_name]]
  current_celltypes <- intersect(current_celltypes, unique(immune_object$major_celltype))

  if (length(current_celltypes) == 0) {
    next
  }

  de_results[[analysis_name]] <- run_seurat_de(
    object = immune_object,
    celltypes = current_celltypes,
    group_label = analysis_name,
    ident_1 = paste(analysis_name, "1+", "pre", sep = "_"),
    ident_2 = paste(analysis_name, "0", "pre", sep = "_"),
    output_prefix = analysis_name
  )

  pseudobulk_results[[analysis_name]] <- run_pseudobulk_de(
    object = immune_object,
    celltypes = current_celltypes,
    output_prefix = analysis_name
  )

}

saveRDS(
  list(
    cell_fraction = cell_fraction,
    paired_fraction_change = paired_fraction_change,
    de_results = de_results,
    pseudobulk_results = pseudobulk_results
  ),
  file = file.path(output_dir, "immune_analysis_results.rds")
)
