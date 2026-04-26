#!/usr/bin/env Rscript

# Reproducible single-cell RNA-seq analysis script for peer review
# Replace the paths and metadata column names below with the corresponding files/columns
# used in the reviewed analysis.

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(cowplot)
  library(readxl)
  library(openxlsx)
})

set.seed(1234)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 200 * 1024^3)

# -------------------------------------------------------------------------
# 1. User-defined inputs
# -------------------------------------------------------------------------

input_object_file <- "data/seurat_object.rds"
clinical_metadata_file <- "data/clinical_metadata.xlsx"
clinical_metadata_sheet <- 1
output_dir <- "peer_review_outputs"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Required metadata columns in clinical_metadata_file:
# sample_id: sample/library identifier matching seurat_obj$orig.ident
# patient_id: patient identifier
# timepoint: sample timepoint, for example "pre" or "post"
# group: comparison group, for example "HER2_0" or "HER2_low"
# treatment: treatment category, for example "IO" or "Chemo"
# tissue: tissue source, for example "Primary" or "LN"
# receptor_status: optional filtering column, for example "TNBC" or "ER_negative"
# repeated_sample: optional indicator where 0 means the sample is retained

# -------------------------------------------------------------------------
# 2. Load data and attach clinical metadata
# -------------------------------------------------------------------------

seurat_obj <- readRDS(input_object_file)

clinical_meta <- read_excel(
  path = clinical_metadata_file,
  sheet = clinical_metadata_sheet
) %>%
  as.data.frame()

stopifnot("sample_id" %in% colnames(clinical_meta))
rownames(clinical_meta) <- clinical_meta$sample_id

metadata_columns <- intersect(
  c(
    "patient_id", "timepoint", "group", "treatment",
    "tissue", "receptor_status", "repeated_sample"
  ),
  colnames(clinical_meta)
)

seurat_obj@meta.data[, metadata_columns] <-
  clinical_meta[seurat_obj@meta.data$orig.ident, metadata_columns, drop = FALSE]

# Keep non-duplicated samples with available group labels when these columns exist.
if (all(c("repeated_sample", "group") %in% colnames(seurat_obj@meta.data))) {
  seurat_obj <- subset(seurat_obj, subset = repeated_sample == 0 & !is.na(group))
}

# -------------------------------------------------------------------------
# 3. Quality control
# -------------------------------------------------------------------------

seurat_obj[["percent_mito"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent_ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
seurat_obj[["percent_hb"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^HB[^(P)]"
)

seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA >= 400 &
    nFeature_RNA <= 8000 &
    nCount_RNA >= 500 &
    percent_mito < 10
)

saveRDS(seurat_obj, file.path(output_dir, "01_qc_filtered_seurat_object.rds"))

# -------------------------------------------------------------------------
# 4. Normalization, integration, dimensionality reduction, and clustering
# -------------------------------------------------------------------------

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

batch_column <- "patient_id"
if (!batch_column %in% colnames(seurat_obj@meta.data)) {
  batch_column <- "orig.ident"
}

seurat_obj <- RunHarmony(seurat_obj, group.by.vars = batch_column)
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:50)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:50)
seurat_obj <- FindClusters(seurat_obj, resolution = 1.0)

saveRDS(seurat_obj, file.path(output_dir, "02_integrated_clustered_seurat_object.rds"))

cluster_plot <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  raster = FALSE
) +
  theme_classic()

ggsave(
  filename = file.path(output_dir, "umap_clusters.pdf"),
  plot = cluster_plot,
  width = 8,
  height = 7
)

# -------------------------------------------------------------------------
# 5. Cell-type annotation
# -------------------------------------------------------------------------
# Edit this table according to marker-gene inspection. The example below is
# intentionally generic and should be replaced with the final cluster labels.

cluster_annotation <- data.frame(
  cluster = levels(seurat_obj$seurat_clusters),
  celltype = levels(seurat_obj$seurat_clusters)
)

# Example:
# cluster_annotation$celltype[cluster_annotation$cluster %in% c("0", "1")] <- "T cell"
# cluster_annotation$celltype[cluster_annotation$cluster %in% c("2")] <- "B cell"
# cluster_annotation$celltype[cluster_annotation$cluster %in% c("3")] <- "Myeloid"
# cluster_annotation$celltype[cluster_annotation$cluster %in% c("4")] <- "Epithelial"

seurat_obj$celltype <- cluster_annotation$celltype[
  match(seurat_obj$seurat_clusters, cluster_annotation$cluster)
]

celltype_plot <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  raster = FALSE
) +
  theme_classic()

ggsave(
  filename = file.path(output_dir, "umap_celltypes.pdf"),
  plot = celltype_plot,
  width = 8,
  height = 7
)

# Optional marker dot plot for final cell-type labels.
marker_list <- list(
  "T cell" = c("CD3D", "CD3E"),
  "B cell" = c("MS4A1", "CD79A"),
  "Plasma cell" = c("MZB1", "JCHAIN"),
  "Myeloid" = c("LYZ", "CD68"),
  "Dendritic cell" = c("FCER1A", "CLEC10A"),
  "Epithelial" = c("EPCAM", "KRT8"),
  "Fibroblast" = c("COL1A1", "DCN"),
  "Endothelial" = c("PECAM1", "VWF")
)

marker_list <- lapply(marker_list, function(x) intersect(x, rownames(seurat_obj)))
marker_list <- marker_list[lengths(marker_list) > 0]

if (length(marker_list) > 0) {
  marker_plot <- DotPlot(seurat_obj, features = marker_list, group.by = "celltype") +
    RotatedAxis() +
    theme_classic() +
    labs(x = NULL, y = NULL)

  ggsave(
    filename = file.path(output_dir, "celltype_marker_dotplot.pdf"),
    plot = marker_plot,
    width = 12,
    height = 5
  )
}

saveRDS(seurat_obj, file.path(output_dir, "03_annotated_seurat_object.rds"))

# -------------------------------------------------------------------------
# 6. Sample-level cell-type fractions
# -------------------------------------------------------------------------

cell_fraction <- seurat_obj@meta.data %>%
  as.data.frame() %>%
  filter(!is.na(celltype)) %>%
  group_by(patient_id, orig.ident, timepoint, group, treatment, tissue, celltype) %>%
  summarise(celltype_n = n(), .groups = "drop") %>%
  group_by(patient_id, orig.ident, timepoint, group, treatment, tissue) %>%
  mutate(total_n = sum(celltype_n), fraction = celltype_n / total_n) %>%
  ungroup()

write.xlsx(
  cell_fraction,
  file = file.path(output_dir, "sample_level_celltype_fractions.xlsx"),
  overwrite = TRUE
)

# Baseline comparison between groups.
baseline_fraction <- cell_fraction %>%
  filter(timepoint == "pre", !is.na(group))

baseline_group_plot <- ggplot(
  baseline_fraction,
  aes(x = celltype, y = fraction, fill = group)
) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  theme_classic() +
  labs(x = NULL, y = "Baseline cell-type fraction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(output_dir, "baseline_celltype_fraction_by_group.pdf"),
  plot = baseline_group_plot,
  width = 12,
  height = 5
)

# Tissue comparison at baseline.
if ("tissue" %in% colnames(baseline_fraction)) {
  tissue_plot <- ggplot(
    baseline_fraction,
    aes(x = celltype, y = fraction, fill = tissue)
  ) +
    geom_boxplot(outlier.shape = NA) +
    stat_compare_means(method = "wilcox.test", label = "p.signif") +
    theme_classic() +
    labs(x = NULL, y = "Baseline cell-type fraction") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(
    filename = file.path(output_dir, "baseline_celltype_fraction_by_tissue.pdf"),
    plot = tissue_plot,
    width = 12,
    height = 5
  )
}

# -------------------------------------------------------------------------
# 7. Paired pre-treatment versus post-treatment comparison
# -------------------------------------------------------------------------

paired_samples <- cell_fraction %>%
  filter(timepoint %in% c("pre", "post")) %>%
  distinct(patient_id, timepoint) %>%
  count(patient_id) %>%
  filter(n == 2) %>%
  pull(patient_id)

paired_fraction <- cell_fraction %>%
  filter(patient_id %in% paired_samples, timepoint %in% c("pre", "post")) %>%
  select(patient_id, group, treatment, timepoint, celltype, fraction) %>%
  tidyr::complete(
    patient_id,
    group,
    treatment,
    timepoint,
    celltype,
    fill = list(fraction = 0)
  )

paired_timepoint_plot <- ggplot(
  paired_fraction,
  aes(x = celltype, y = fraction, fill = timepoint)
) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test", paired = TRUE, label = "p.signif") +
  theme_classic() +
  labs(x = NULL, y = "Cell-type fraction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(output_dir, "paired_pre_post_celltype_fraction.pdf"),
  plot = paired_timepoint_plot,
  width = 12,
  height = 5
)

fraction_change <- paired_fraction %>%
  select(patient_id, group, treatment, celltype, timepoint, fraction) %>%
  pivot_wider(names_from = timepoint, values_from = fraction, values_fill = 0) %>%
  mutate(fraction_change = post - pre)

write.xlsx(
  fraction_change,
  file = file.path(output_dir, "paired_celltype_fraction_change.xlsx"),
  overwrite = TRUE
)

change_plot <- ggplot(
  fraction_change,
  aes(x = celltype, y = fraction_change, fill = group)
) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  theme_classic() +
  labs(x = NULL, y = "Post-treatment minus pre-treatment fraction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(output_dir, "paired_celltype_fraction_change_by_group.pdf"),
  plot = change_plot,
  width = 12,
  height = 5
)

# -------------------------------------------------------------------------
# 8. Session information
# -------------------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "session_info.txt")
)
