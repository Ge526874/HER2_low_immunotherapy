############################################################
## Clinical analysis script for peer review
## Purpose: HER2 subgroup comparison, pCR analysis, and PSM sensitivity analyses
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(tableone)
  library(MatchIt)
})

## ---- 1. Load clinical data ----
## Replace the input file with the reviewer-accessible clinical data file.
clinical_data <- read.csv("clinical_data.csv", check.names = FALSE)

## ---- 2. Standardize variable names ----
## Please adapt the column names on the right if the shared dataset uses different names.
analysis_data <- clinical_data %>%
  transmute(
    her2_ihc = HER2_IHC,
    pcr = pCR,
    treatment_group = treatment_group,
    immunotherapy_type = immunotherapy_type,
    age = age,
    grade = grade,
    ki67 = Ki67,
    pdl1_cps = PDL1_CPS,
    clinical_t = clinical_T,
    clinical_n = clinical_N,
    efs_event = EFS_event,
    efs_time = EFS_time
  ) %>%
  mutate(
    her2_group = if_else(her2_ihc == "0", "HER2-0", "HER2-low"),
    her2_binary = if_else(her2_ihc == "0", 0, 1),
    pdl1_available = if_else(is.na(pdl1_cps), 0, 1),
    immunotherapy = treatment_group == "immunotherapy",
    pcr = factor(pcr),
    her2_group = factor(her2_group, levels = c("HER2-0", "HER2-low")),
    clinical_t = factor(clinical_t),
    clinical_n = factor(clinical_n)
  )

## ---- 3. Utility functions ----
summarize_pcr_by_her2 <- function(data) {
  summary_table <- data %>%
    filter(!is.na(her2_group), !is.na(pcr)) %>%
    count(her2_group, pcr, name = "n") %>%
    group_by(her2_group) %>%
    mutate(
      percent = 100 * n / sum(n),
      label = paste0(round(percent, 1), "%")
    ) %>%
    ungroup()

  contingency_table <- table(data$her2_group, data$pcr)
  fisher_result <- fisher.test(contingency_table)

  list(
    summary_table = summary_table,
    contingency_table = contingency_table,
    fisher_test = fisher_result
  )
}

plot_pcr_by_her2 <- function(summary_table, p_value, output_file) {
  p <- ggplot(summary_table, aes(x = her2_group, y = percent, fill = pcr)) +
    geom_col(width = 0.7) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 4
    ) +
    annotate(
      "text",
      x = 1.5,
      y = 105,
      label = paste0("Fisher's exact test, P = ", format.pval(p_value, digits = 3, eps = 0.001)),
      size = 4
    ) +
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
    labs(x = NULL, y = "Percentage", fill = "pCR") +
    theme_classic(base_size = 12)

  ggsave(output_file, p, width = 4, height = 4)
  p
}

impute_psm_covariates <- function(data) {
  data %>%
    mutate(
      pdl1_cps = if_else(is.na(pdl1_cps), median(pdl1_cps, na.rm = TRUE), pdl1_cps)
    )
}

perform_psm <- function(data) {
  psm_formula <- her2_binary ~ grade + ki67 + pdl1_cps + clinical_t + clinical_n + age
  }

  matchit(
    formula = psm_formula,
    data = data,
    method = "nearest",
    caliper = 0.2,
    ratio = 1
  )
}

analyze_matched_pcr <- function(matched_data, output_file) {
  pcr_result <- summarize_pcr_by_her2(matched_data)
  plot_pcr_by_her2(
    summary_table = pcr_result$summary_table,
    p_value = pcr_result$fisher_test$p.value,
    output_file = output_file
  )

  pcr_result
}

## ---- 4. Descriptive statistics ----
vars_for_table <- c(
  "age", "grade", "ki67", "pdl1_cps",
  "clinical_t", "clinical_n", "pcr", "efs_event", "efs_time"
)

table_immunotherapy <- CreateTableOne(
  vars = vars_for_table,
  strata = "her2_group",
  data = filter(analysis_data, immunotherapy),
  test = TRUE
)

table_non_immunotherapy <- CreateTableOne(
  vars = vars_for_table,
  strata = "her2_group",
  data = filter(analysis_data, !immunotherapy),
  test = TRUE
)

print(table_immunotherapy, showAllLevels = TRUE)
print(table_non_immunotherapy, showAllLevels = TRUE)

## ---- 5. Unmatched pCR comparison ----
unmatched_immunotherapy <- summarize_pcr_by_her2(
  filter(analysis_data, immunotherapy)
)

unmatched_non_immunotherapy <- summarize_pcr_by_her2(
  filter(analysis_data, !immunotherapy)
)

plot_pcr_by_her2(
  unmatched_immunotherapy$summary_table,
  unmatched_immunotherapy$fisher_test$p.value,
  "unmatched_pcr_immunotherapy.pdf"
)

plot_pcr_by_her2(
  unmatched_non_immunotherapy$summary_table,
  unmatched_non_immunotherapy$fisher_test$p.value,
  "unmatched_pcr_non_immunotherapy.pdf"
)

## ---- 6. Propensity score matching ----
psm_data <- analysis_data %>%
  filter(!is.na(her2_binary), !is.na(pcr)) %>%
  impute_psm_covariates()

psm_immunotherapy <- perform_psm(
  filter(psm_data, immunotherapy)
)

psm_non_immunotherapy <- perform_psm(
  filter(psm_data, !immunotherapy)
)

matched_immunotherapy <- match.data(psm_immunotherapy)
matched_non_immunotherapy <- match.data(psm_non_immunotherapy)

matched_result_immunotherapy <- analyze_matched_pcr(
  matched_immunotherapy,
  "matched_pcr_immunotherapy.pdf"
)

matched_result_non_immunotherapy <- analyze_matched_pcr(
  matched_non_immunotherapy,
  "matched_pcr_non_immunotherapy.pdf"
)

## ---- 7. Export key result tables ----
write.csv(
  unmatched_immunotherapy$summary_table,
  "unmatched_pcr_immunotherapy_summary.csv",
  row.names = FALSE
)

write.csv(
  unmatched_non_immunotherapy$summary_table,
  "unmatched_pcr_non_immunotherapy_summary.csv",
  row.names = FALSE
)

write.csv(
  matched_result_immunotherapy$summary_table,
  "matched_pcr_immunotherapy_summary.csv",
  row.names = FALSE
)

write.csv(
  matched_result_non_immunotherapy$summary_table,
  "matched_pcr_non_immunotherapy_summary.csv",
  row.names = FALSE
)

sessionInfo()
