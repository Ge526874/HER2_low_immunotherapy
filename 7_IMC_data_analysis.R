# Peer-review R script for spatial single-cell analysis

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
  library(corrplot)
  library(ggsci)
  library(sf)
})

set.seed(123)

input_cell_file <- "cell_metadata.csv"
input_clinical_file <- "clinical_metadata.csv"
input_neighbor_file <- "neighbor_interactions.csv"
output_dir <- "Figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cell_data <- read.csv(input_cell_file, check.names = FALSE)
clinical_data <- read.csv(input_clinical_file, check.names = FALSE)
neighbor_data <- read.csv(input_neighbor_file, check.names = FALSE)

clinical_data <- clinical_data %>%
  dplyr::select(patient_id, treatment_arm, erbb2_score, per_protocol, response)

cell_data <- cell_data %>%
  left_join(clinical_data, by = "patient_id")

marker_annotation <- data.frame(
  marker = c(
    "CD11c", "CD15", "CD163", "CD68", "MPO",
    "Caveolin.1", "CD31", "PDPN", "PDGFRB", "SMA", "Calponin", "Vimentin",
    "PD.L1..SP142.", "PD.L1..73.10.", "IDO", "PD.1", "OX40", "ICOS",
    "CK5.14", "CK8.18", "panCK", "AR", "CA9", "c.PARP", "pH2AX", "Ki67",
    "CD20", "CD79a", "CD56", "CD3", "CD4", "CD8", "FOXP3", "GATA3", "Helios",
    "T.bet", "TCF1", "TOX", "GZMB", "CD45", "HLA.ABC", "HLA.DR"
  ),
  marker_group = c(
    rep("Myeloid", 5), rep("Mesenchymal", 7), rep("Checkpoint", 6),
    rep("Epithelial", 4), "Hypoxia", rep("Cell_state", 3),
    rep("Lymphoid", 13), "Pan_immune", rep("MHC", 2)
  ),
  tumor_marker = c(
    0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1,
    1, 1, 1, 0, 0, 0,
    1, 1, 1, 1, 1, 1, 1, 1,
    0, 0, 1, 0, 0, 0, 0, 1, 1,
    0, 1, 0, 1, 0, 1, 1
  )
)

cell_type_annotation <- data.frame(
  cell_label = c(
    "Myofibroblasts", "Fibroblasts", "Endothelial", "PDPN^+Stromal",
    "CD4^+TCF1^+T", "Treg", "CD4^+PD1^+T",
    "CD8^+T", "CD8^+TCF1^+T", "CD8^+PD1^+T_{Ex}", "CD8^+GZMB^+T",
    "CD56^+NK", "CD20^+B", "CD79a^+Plasma",
    "DCs", "M2 Mac", "PD-L1^+IDO^+APCs", "PD-L1^+APCs", "Neutrophils", "CA9^+"
  ),
  major_cell_type = c(
    rep("Stromal cells", 4), rep("CD4+ T cells", 3), rep("CD8+ T cells", 4),
    "NK cells", rep("B cells", 2), rep("APC cells", 4), rep("Other cells", 2)
  )
)

safe_cor_test <- function(x, y) {
  valid <- complete.cases(x, y)
  if (sum(valid) < 3 || length(unique(x[valid])) < 2 || length(unique(y[valid])) < 2) {
    return(data.frame(correlation = NA_real_, p_value = NA_real_))
  }
  test_result <- suppressWarnings(cor.test(x[valid], y[valid], method = "spearman"))
  data.frame(correlation = unname(test_result$estimate), p_value = test_result$p.value)
}

plot_marker_correlation <- function(plot_data, x_var, y_var, facet_var, x_label, y_label, output_file, width, height) {
  p <- ggplot(plot_data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(size = 1.5) +
    geom_smooth(method = "lm", se = TRUE) +
    facet_wrap(stats::as.formula(paste("~", facet_var)), scales = "free") +
    stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~~"))) +
    labs(x = x_label, y = y_label) +
    theme_classic()
  ggsave(file.path(output_dir, output_file), p, width = width, height = height)
  p
}

plot_correlation_heatmap <- function(correlation_matrix, pvalue_matrix, output_file, width = 6, height = 15) {
  pdf(file.path(output_dir, output_file), width = width, height = height)
  corrplot(
    correlation_matrix,
    method = "square",
    p.mat = pvalue_matrix,
    insig = "label_sig",
    pch.cex = 1.5,
    col = rev(COL2("PuOr", 200)),
    sig.level = c(0.001, 0.01, 0.05),
    tl.col = "black",
    na.label = " "
  )
  dev.off()
}

tumor_markers <- marker_annotation %>%
  filter(tumor_marker == 1) %>%
  pull(marker)

tumor_baseline <- cell_data %>%
  filter(
    biopsy_phase == "Baseline",
    !is.na(erbb2_score),
    compartment %in% c("DCIS", "invasive", "normal")
  ) %>%
  mutate(erbb2_score_for_analysis = erbb2_score)

tumor_marker_long <- tumor_baseline %>%
  dplyr::select(patient_id, erbb2_score_for_analysis, all_of(tumor_markers)) %>%
  pivot_longer(cols = all_of(tumor_markers), names_to = "marker", values_to = "marker_level") %>%
  left_join(marker_annotation, by = "marker")

tumor_patient_marker <- tumor_marker_long %>%
  group_by(patient_id, erbb2_score_for_analysis, marker, marker_group) %>%
  summarise(mean_marker_level = mean(marker_level, na.rm = TRUE), .groups = "drop")

selected_tumor_markers <- c("PD.L1..73.10.", "PD.L1..SP142.", "IDO", "HLA.ABC", "HLA.DR", "Vimentin", "CD56", "AR", "CA9")
plot_marker_correlation(
  tumor_patient_marker %>% filter(marker %in% selected_tumor_markers),
  x_var = "erbb2_score_for_analysis",
  y_var = "mean_marker_level",
  facet_var = "marker",
  x_label = "Predicted ERBB2 score",
  y_label = "Mean marker level",
  output_file = "tumor_marker_baseline_correlation.pdf",
  width = 9,
  height = 9
)

paired_patients <- cell_data %>%
  filter(biopsy_phase %in% c("Baseline", "On-treatment"), !is.na(erbb2_score)) %>%
  distinct(patient_id, biopsy_phase) %>%
  count(patient_id) %>%
  filter(n == 2) %>%
  pull(patient_id)

tumor_dynamic <- cell_data %>%
  filter(
    biopsy_phase %in% c("Baseline", "On-treatment"),
    patient_id %in% paired_patients,
    !is.na(erbb2_score),
    compartment %in% c("DCIS", "invasive", "normal")
  ) %>%
  mutate(erbb2_score_for_analysis = erbb2_score)

tumor_dynamic_long <- tumor_dynamic %>%
  dplyr::select(patient_id, erbb2_score_for_analysis, biopsy_phase, treatment_arm, per_protocol, all_of(tumor_markers)) %>%
  pivot_longer(cols = all_of(tumor_markers), names_to = "marker", values_to = "marker_level") %>%
  left_join(marker_annotation, by = "marker")

tumor_dynamic_change <- tumor_dynamic_long %>%
  group_by(patient_id, erbb2_score_for_analysis, biopsy_phase, treatment_arm, per_protocol, marker, marker_group) %>%
  summarise(mean_marker_level = mean(marker_level, na.rm = TRUE), .groups = "drop") %>%
  filter(per_protocol == TRUE) %>%
  dplyr::select(patient_id, erbb2_score_for_analysis, biopsy_phase, treatment_arm, marker, marker_group, mean_marker_level) %>%
  pivot_wider(names_from = biopsy_phase, values_from = mean_marker_level, values_fill = 0) %>%
  mutate(marker_level_change = `On-treatment` - Baseline)

plot_marker_correlation(
  tumor_dynamic_change,
  x_var = "erbb2_score_for_analysis",
  y_var = "marker_level_change",
  facet_var = "marker",
  x_label = "Predicted ERBB2 score",
  y_label = "On-treatment minus baseline marker level",
  output_file = "tumor_marker_dynamic_change.pdf",
  width = 15,
  height = 12
)

tme_markers <- marker_annotation %>%
  filter(marker_group != "Epithelial") %>%
  pull(marker)

tme_baseline <- cell_data %>%
  filter(biopsy_phase == "Baseline", !is.na(erbb2_score), compartment == "TME") %>%
  left_join(cell_type_annotation, by = "cell_label") %>%
  mutate(erbb2_score_for_analysis = erbb2_score)

tme_marker_long <- tme_baseline %>%
  dplyr::select(patient_id, erbb2_score_for_analysis, major_cell_type, all_of(tme_markers)) %>%
  pivot_longer(cols = all_of(tme_markers), names_to = "marker", values_to = "marker_level") %>%
  left_join(marker_annotation, by = "marker")

tme_patient_marker <- tme_marker_long %>%
  group_by(patient_id, erbb2_score_for_analysis, major_cell_type, marker, marker_group) %>%
  summarise(mean_marker_level = mean(marker_level, na.rm = TRUE), .groups = "drop")

tme_cor_results <- tme_patient_marker %>%
  group_by(major_cell_type, marker, marker_group) %>%
  group_modify(~ safe_cor_test(.x$erbb2_score_for_analysis, .x$mean_marker_level)) %>%
  ungroup()

tme_cor_matrix <- tme_cor_results %>%
  dplyr::select(marker, major_cell_type, correlation) %>%
  pivot_wider(names_from = marker, values_from = correlation) %>%
  column_to_rownames("major_cell_type") %>%
  as.matrix()

tme_p_matrix <- tme_cor_results %>%
  dplyr::select(marker, major_cell_type, p_value) %>%
  pivot_wider(names_from = marker, values_from = p_value) %>%
  column_to_rownames("major_cell_type") %>%
  as.matrix()

plot_correlation_heatmap(t(tme_cor_matrix), t(tme_p_matrix), "tme_marker_correlation_heatmap.pdf")

cell_fraction <- cell_data %>%
  filter(!is.na(erbb2_score), !(compartment %in% c("DCIS", "invasive", "normal"))) %>%
  group_by(patient_id, biopsy_phase) %>%
  mutate(total_cells = n()) %>%
  ungroup() %>%
  group_by(patient_id, biopsy_phase, erbb2_score, per_protocol, compartment, cell_label) %>%
  summarise(cell_fraction = n() / unique(total_cells), .groups = "drop") %>%
  mutate(erbb2_score_for_analysis = erbb2_score)

plot_marker_correlation(
  cell_fraction %>% filter(biopsy_phase == "Baseline"),
  x_var = "erbb2_score_for_analysis",
  y_var = "cell_fraction",
  facet_var = "cell_label",
  x_label = "Predicted ERBB2 score",
  y_label = "Cell fraction",
  output_file = "baseline_cell_fraction_correlation.pdf",
  width = 15,
  height = 15
)

cell_fraction_change <- cell_fraction %>%
  filter(patient_id %in% paired_patients, per_protocol == TRUE) %>%
  dplyr::select(patient_id, erbb2_score_for_analysis, biopsy_phase, cell_label, cell_fraction) %>%
  pivot_wider(names_from = biopsy_phase, values_from = cell_fraction, values_fill = 0) %>%
  mutate(fraction_change = `On-treatment` - Baseline)

plot_marker_correlation(
  cell_fraction_change,
  x_var = "erbb2_score_for_analysis",
  y_var = "fraction_change",
  facet_var = "cell_label",
  x_label = "Predicted ERBB2 score",
  y_label = "On-treatment minus baseline cell fraction",
  output_file = "cell_fraction_dynamic_change.pdf",
  width = 15,
  height = 9
)

ki67_by_cell_type <- cell_data %>%
  left_join(cell_type_annotation, by = "cell_label") %>%
  filter(!is.na(erbb2_score), !is.na(major_cell_type)) %>%
  group_by(patient_id, biopsy_phase, erbb2_score, major_cell_type) %>%
  summarise(mean_ki67_level = mean(Ki67, na.rm = TRUE), .groups = "drop") %>%
  mutate(erbb2_score_for_analysis = erbb2_score)

plot_marker_correlation(
  ki67_by_cell_type %>% filter(biopsy_phase == "Baseline", major_cell_type == "CD8+ T cells"),
  x_var = "erbb2_score_for_analysis",
  y_var = "mean_ki67_level",
  facet_var = "major_cell_type",
  x_label = "Predicted ERBB2 score",
  y_label = "Mean Ki67 level",
  output_file = "cd8_t_cell_ki67_correlation.pdf",
  width = 3.5,
  height = 3
)

get_tissue_area <- function(x, y) {
  coordinate_data <- data.frame(x = x, y = y) %>% distinct()
  if (nrow(coordinate_data) < 3) return(NA_real_)
  hull_indices <- chull(coordinate_data$x, coordinate_data$y)
  hull_indices <- c(hull_indices, hull_indices[1])
  hull_points <- coordinate_data[hull_indices, ]
  as.numeric(st_area(st_polygon(list(as.matrix(hull_points)))))
}

patient_tissue_area <- cell_data %>%
  group_by(patient_id, biopsy_phase, image_id) %>%
  summarise(tissue_area = get_tissue_area(x_coordinate, y_coordinate), .groups = "drop") %>%
  group_by(patient_id, biopsy_phase) %>%
  summarise(patient_tissue_area = sum(tissue_area, na.rm = TRUE), .groups = "drop")

cell_data_with_area <- cell_data %>%
  left_join(patient_tissue_area, by = c("patient_id", "biopsy_phase"))

cell_density <- cell_data_with_area %>%
  filter(!is.na(erbb2_score)) %>%
  group_by(patient_id, biopsy_phase, erbb2_score, cell_label) %>%
  summarise(cell_density = n() / unique(patient_tissue_area), .groups = "drop") %>%
  mutate(erbb2_score_for_analysis = erbb2_score)

plot_marker_correlation(
  cell_density %>% filter(biopsy_phase == "Baseline"),
  x_var = "erbb2_score_for_analysis",
  y_var = "cell_density",
  facet_var = "cell_label",
  x_label = "Predicted ERBB2 score",
  y_label = "Cell density",
  output_file = "baseline_cell_density_correlation.pdf",
  width = 15,
  height = 15
)

image_cell_counts <- cell_data %>%
  group_by(patient_id, biopsy_phase, image_id) %>%
  summarise(image_cell_count = n(), .groups = "drop")

patient_cell_counts <- cell_data %>%
  group_by(patient_id, biopsy_phase) %>%
  summarise(patient_cell_count = n(), .groups = "drop")

epithelial_neighbor_counts <- neighbor_data %>%
  filter(from_epithelial == TRUE) %>%
  group_by(image_id, neighbor_cell_label) %>%
  summarise(interactions = n(), .groups = "drop") %>%
  left_join(image_cell_counts, by = "image_id") %>%
  left_join(patient_cell_counts, by = c("patient_id", "biopsy_phase")) %>%
  group_by(patient_id, biopsy_phase, neighbor_cell_label) %>%
  summarise(interaction_index = sum(interactions, na.rm = TRUE) / unique(patient_cell_count), .groups = "drop") %>%
  left_join(clinical_data, by = "patient_id") %>%
  filter(!is.na(erbb2_score)) %>%
  mutate(erbb2_score_for_analysis = log2(pmax(erbb2_score, 0) + 1))

immune_neighbor_labels <- c(
  "M2 Mac", "DCs", "Neutrophils", "PD-L1^+IDO^+APCs", "PD-L1^+APCs",
  "CD4^+TCF1^+T", "CD20^+B", "CD8^+T", "CD79a^+Plasma", "CD8^+TCF1^+T",
  "CD8^+PD1^+T_{Ex}", "CD8^+GZMB^+T", "Treg", "CD4^+PD1^+T", "CD56^+NK"
)

plot_marker_correlation(
  epithelial_neighbor_counts %>% filter(biopsy_phase == "Baseline", neighbor_cell_label %in% immune_neighbor_labels),
  x_var = "erbb2_score_for_analysis",
  y_var = "interaction_index",
  facet_var = "neighbor_cell_label",
  x_label = "Log2-transformed predicted ERBB2 score",
  y_label = "Epithelial-neighbor interaction index",
  output_file = "epithelial_immune_neighbor_correlation.pdf",
  width = 18,
  height = 11
)

write.csv(tumor_patient_marker, file.path(output_dir, "tumor_marker_patient_summary.csv"), row.names = FALSE)
write.csv(tumor_dynamic_change, file.path(output_dir, "tumor_marker_dynamic_change.csv"), row.names = FALSE)
write.csv(tme_cor_results, file.path(output_dir, "tme_marker_correlation_results.csv"), row.names = FALSE)
write.csv(cell_fraction, file.path(output_dir, "cell_fraction_summary.csv"), row.names = FALSE)
write.csv(cell_density, file.path(output_dir, "cell_density_summary.csv"), row.names = FALSE)
write.csv(epithelial_neighbor_counts, file.path(output_dir, "epithelial_neighbor_interaction_summary.csv"), row.names = FALSE)
