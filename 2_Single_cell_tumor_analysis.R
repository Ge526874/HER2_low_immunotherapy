options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(readxl)
  library(DESeq2)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
})

set.seed(1234)

# Input files. Replace these paths with reviewer-accessible files.
seurat_object_file <- "input/annotated_scRNAseq_seurat_object.rds"
clinical_metadata_file <- "input/sample_clinical_metadata.xlsx"
output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Required clinical metadata columns:
# sample_id, timepoint, treatment_group, HER2_status, ER_status, pathological_response, batch
sc_object <- readRDS(seurat_object_file)
clinical_metadata <- read_excel(clinical_metadata_file) |> as.data.frame()
rownames(clinical_metadata) <- clinical_metadata$sample_id

sc_object$timepoint <- clinical_metadata[sc_object$orig.ident, "timepoint"]
sc_object$treatment_group <- clinical_metadata[sc_object$orig.ident, "treatment_group"]
sc_object$HER2_status <- clinical_metadata[sc_object$orig.ident, "HER2_status"]
sc_object$ER_status <- clinical_metadata[sc_object$orig.ident, "ER_status"]
sc_object$pathological_response <- clinical_metadata[sc_object$orig.ident, "pathological_response"]
sc_object$batch <- clinical_metadata[sc_object$orig.ident, "batch"]

Idents(sc_object) <- "celltype"

# UMAP visualization of major cell types.
p_umap <- DimPlot(sc_object, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) +
  theme_classic() +
  labs(title = NULL)
ggsave(file.path(output_dir, "major_celltype_umap.pdf"), p_umap, width = 7, height = 6)

# Patient-level ERBB2 expression in epithelial/tumor cells at baseline.
DefaultAssay(sc_object) <- "RNA"
expression_matrix <- GetAssayData(sc_object, assay = "RNA", slot = "data")
sc_object$ERBB2_expression <- as.numeric(expression_matrix["ERBB2", colnames(sc_object)])

epithelial_object <- subset(sc_object, subset = celltype %in% c("Epithelial", "Tumor"))
epithelial_metadata <- epithelial_object@meta.data

patient_erbb2 <- epithelial_metadata |>
  filter(timepoint == "pre") |>
  group_by(orig.ident, HER2_status) |>
  summarise(mean_ERBB2 = mean(ERBB2_expression, na.rm = TRUE), .groups = "drop")

p_erbb2 <- ggplot(patient_erbb2, aes(x = HER2_status, y = mean_ERBB2, fill = HER2_status)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.12, size = 1.8) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_classic() +
  labs(x = "HER2 status", y = "Mean ERBB2 expression") +
  theme(legend.position = "none")
ggsave(file.path(output_dir, "baseline_patient_level_ERBB2.pdf"), p_erbb2, width = 3.2, height = 4)

# Patient-level expression of immune-related genes in epithelial/tumor cells.
genes_of_interest <- c("ERBB2", "HLA-A", "HLA-B", "HLA-C", "B2M", "CD274", "PDCD1LG2")
genes_of_interest <- intersect(genes_of_interest, rownames(expression_matrix))

for (gene in genes_of_interest) {
  epithelial_metadata[[paste0(gene, "_expression")]] <- as.numeric(GetAssayData(epithelial_object, assay = "RNA", slot = "data")[gene, colnames(epithelial_object)])
}

patient_gene_expression <- epithelial_metadata |>
  filter(timepoint == "pre") |>
  group_by(orig.ident, HER2_status) |>
  summarise(across(ends_with("_expression"), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  pivot_longer(cols = ends_with("_expression"), names_to = "gene", values_to = "mean_expression") |>
  mutate(gene = sub("_expression$", "", gene))

p_gene_expression <- ggplot(patient_gene_expression, aes(x = HER2_status, y = mean_expression, fill = HER2_status)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.12, size = 0.8) +
  facet_wrap(~ gene, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_classic() +
  labs(x = "HER2 status", y = "Mean expression") +
  theme(legend.position = "none")
ggsave(file.path(output_dir, "baseline_patient_level_marker_genes.pdf"), p_gene_expression, width = 8, height = 5)

# Differential expression between HER2 groups within epithelial/tumor cells.
epithelial_object_baseline <- subset(epithelial_object, subset = timepoint == "pre" & !is.na(HER2_status))
Idents(epithelial_object_baseline) <- "HER2_status"

de_results <- FindMarkers(
  epithelial_object_baseline,
  ident.1 = "HER2-low",
  ident.2 = "HER2-0",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.1
)
de_results$gene <- rownames(de_results)
write.csv(de_results, file.path(output_dir, "baseline_epithelial_DE_HER2_low_vs_zero.csv"), row.names = FALSE)

# Pseudobulk differential expression at the sample level.
pseudobulk_object <- AggregateExpression(
  epithelial_object_baseline,
  assays = "RNA",
  return.seurat = TRUE,
  group.by = c("orig.ident", "HER2_status")
)

pseudobulk_counts <- GetAssayData(pseudobulk_object, assay = "RNA", slot = "counts")
pseudobulk_metadata <- pseudobulk_object@meta.data
pseudobulk_metadata$sample_id <- sub("_.*$", "", rownames(pseudobulk_metadata))
pseudobulk_metadata$HER2_status <- factor(pseudobulk_metadata$HER2_status)

common_samples <- intersect(colnames(pseudobulk_counts), rownames(pseudobulk_metadata))
pseudobulk_counts <- pseudobulk_counts[, common_samples]
pseudobulk_metadata <- pseudobulk_metadata[common_samples, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = round(pseudobulk_counts),
  colData = pseudobulk_metadata,
  design = ~ HER2_status
)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
pseudobulk_de <- results(dds, contrast = c("HER2_status", "HER2-low", "HER2-0")) |> as.data.frame()
pseudobulk_de$gene <- rownames(pseudobulk_de)
write.csv(pseudobulk_de, file.path(output_dir, "baseline_epithelial_pseudobulk_DE_HER2_low_vs_zero.csv"), row.names = FALSE)

saveRDS(sc_object, file.path(output_dir, "analysis_ready_scRNAseq_object.rds"))
