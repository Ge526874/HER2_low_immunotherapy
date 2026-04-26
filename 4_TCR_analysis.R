options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(scRepertoire)
  library(SingleCellExperiment)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
})

set.seed(123)

# Input files expected from the user:
# 1) A folder containing one subfolder per sample, each with all_contig_annotations.csv.
# 2) A sample metadata table with sample_id, patient_id, timepoint, treatment_group, marker_group, and optional response/status columns.
# 3) Optional: a SingleCellExperiment object with cell-level TCR clone annotation for cell-type-specific analyses.

contig_parent_dir <- "path/to/tcr_contig_folders"
metadata_file <- "path/to/sample_metadata.csv"
sce_file <- "path/to/single_cell_object.rds"
output_dir <- "peer_review_tcr_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

sample_metadata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)
sample_metadata$sample_id <- rownames(sample_metadata)

sample_metadata <- sample_metadata %>%
  mutate(
    timepoint = factor(timepoint, levels = c("pre", "post")),
    group_label = paste(marker_group, treatment_group, sep = " | "),
    sample_group = paste(marker_group, treatment_group, timepoint, sep = "_")
  )

sample_ids <- sample_metadata$sample_id
contig_list <- lapply(sample_ids, function(sample_id) {
  read.csv(file.path(contig_parent_dir, sample_id, "all_contig_annotations.csv"), check.names = FALSE)
})
names(contig_list) <- sample_ids

combined_tcr <- combineTCR(
  contig_list,
  samples = sample_ids,
  removeNA = FALSE,
  removeMulti = FALSE,
  filterMulti = FALSE
)

combined_tcr <- addVariable(combined_tcr, "timepoint", sample_metadata$timepoint)
combined_tcr <- addVariable(combined_tcr, "treatment_group", sample_metadata$treatment_group)
combined_tcr <- addVariable(combined_tcr, "marker_group", sample_metadata$marker_group)
combined_tcr <- addVariable(combined_tcr, "patient_id", sample_metadata$patient_id)
combined_tcr <- addVariable(combined_tcr, "group_label", sample_metadata$group_label)
combined_tcr <- addVariable(combined_tcr, "sample_group", sample_metadata$sample_group)

saveRDS(combined_tcr, file.path(output_dir, "combined_tcr.rds"))

basic_clonal_quantification <- clonalQuant(
  combined_tcr,
  cloneCall = "aa",
  chain = "both",
  scale = TRUE
)
ggsave(file.path(output_dir, "clonal_quantification.pdf"), basic_clonal_quantification, width = 7, height = 5)

clonal_homeostasis_plot <- clonalHomeostasis(
  combined_tcr,
  cloneSize = c(Rare = 3e-05, Small = 1e-04, Medium = 0.001, Large = 0.01, Hyperexpanded = 1),
  group.by = "sample_group",
  cloneCall = "aa"
) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(output_dir, "clonal_homeostasis_by_group.pdf"), clonal_homeostasis_plot, width = 10, height = 6)

shannon_diversity <- clonalDiversity(
  combined_tcr,
  cloneCall = "gene",
  group.by = "sample",
  n.boots = 100
)

shannon_df <- shannon_diversity$data %>%
  left_join(sample_metadata, by = c("sample" = "sample_id")) %>%
  filter(timepoint %in% c("pre", "post")) %>%
  group_by(group_label, patient_id) %>%
  filter(n_distinct(timepoint) == 2) %>%
  ungroup()

write.csv(shannon_df, file.path(output_dir, "shannon_diversity_by_sample.csv"), row.names = FALSE)

shannon_plot <- ggplot(shannon_df, aes(x = timepoint, y = value, fill = timepoint)) +
  facet_wrap(~ group_label, scales = "free", ncol = 4) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  geom_line(aes(group = patient_id), alpha = 0.3) +
  stat_compare_means(method = "wilcox.test", paired = TRUE, label = "p.format") +
  theme_classic() +
  labs(x = NULL, y = "Shannon diversity") +
  theme(strip.background = element_blank())
ggsave(file.path(output_dir, "paired_shannon_diversity.pdf"), shannon_plot, width = 10, height = 4)

if (file.exists(sce_file)) {
  sce_object <- readRDS(sce_file)
  cell_metadata <- as.data.frame(colData(sce_object))

  clone_frequency <- cell_metadata %>%
    filter(!is.na(sample_id), !is.na(cell_type), !is.na(tcr_clone_id)) %>%
    group_by(sample_id, cell_type, tcr_clone_id) %>%
    summarise(cell_count = n(), .groups = "drop") %>%
    group_by(sample_id, cell_type) %>%
    mutate(clone_frequency = cell_count / sum(cell_count)) %>%
    ungroup() %>%
    left_join(sample_metadata, by = "sample_id") %>%
    filter(timepoint %in% c("pre", "post")) %>%
    group_by(group_label, patient_id) %>%
    filter(n_distinct(timepoint) == 2) %>%
    ungroup()

  clone_frequency_wide <- clone_frequency %>%
    select(patient_id, cell_type, marker_group, treatment_group, timepoint, tcr_clone_id, clone_frequency) %>%
    pivot_wider(
      names_from = timepoint,
      values_from = clone_frequency,
      values_fill = 0
    ) %>%
    mutate(
      distance_to_diagonal = abs(post - pre) / sqrt(2),
      signed_distance = (post - pre) / sqrt(2),
      clone_weight = pre + post
    )

  expansion_score <- clone_frequency_wide %>%
    group_by(patient_id, cell_type, marker_group, treatment_group) %>%
    summarise(
      mean_distance = mean(distance_to_diagonal),
      median_distance = median(distance_to_diagonal),
      weighted_distance = sum(distance_to_diagonal * clone_weight) / sum(clone_weight),
      mean_signed_distance = mean(signed_distance),
      weighted_signed_distance = sum(signed_distance * clone_weight) / sum(clone_weight),
      .groups = "drop"
    ) %>%
    mutate(group_label = paste(marker_group, treatment_group, sep = " | "))

  write.csv(expansion_score, file.path(output_dir, "cell_type_clone_expansion_score.csv"), row.names = FALSE)

  expansion_plot <- ggplot(expansion_score, aes(x = group_label, y = median_distance, fill = group_label)) +
    facet_wrap(~ cell_type, scales = "free", ncol = 4) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.6) +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    theme_classic() +
    labs(x = NULL, y = "Median distance to diagonal") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")
  ggsave(file.path(output_dir, "cell_type_clone_expansion_score.pdf"), expansion_plot, width = 12, height = 10)
}
