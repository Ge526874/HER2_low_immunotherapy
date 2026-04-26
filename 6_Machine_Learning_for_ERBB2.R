# Peer-review R script: ERBB2 expression prediction with cross-validation and hyperparameter selection

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readr)
  library(mlr3)
  library(mlr3learners)
  library(mlr3pipelines)
  library(mlr3extralearners)
  library(glmnet)
  library(ranger)
  library(ggplot2)
  library(ggpubr)
})

set.seed(123)
lgr::get_logger("mlr3")$set_threshold("warn")

output_dir <- "peer_review_outputs"
model_dir <- file.path(output_dir, "models")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

target_gene <- "ERBB2"

read_expression_matrix <- function(file_path) {
  expression_data <- read_csv(file_path, show_col_types = FALSE)
  expression_data <- as.data.frame(expression_data)
  rownames(expression_data) <- expression_data[[1]]
  expression_data <- expression_data[, -1, drop = FALSE]
  return(expression_data)
}

prepare_model_matrix <- function(expression_matrix, feature_genes, target_gene = "ERBB2", log_transform = TRUE) {
  expression_matrix <- expression_matrix[!duplicated(rownames(expression_matrix)), , drop = FALSE]
  expression_matrix <- expression_matrix[, colSums(is.na(expression_matrix)) < nrow(expression_matrix), drop = FALSE]

  selected_genes <- intersect(c(feature_genes, target_gene), rownames(expression_matrix))
  model_data <- t(expression_matrix[selected_genes, , drop = FALSE])
  model_data <- as.data.frame(model_data)

  if (log_transform) {
    model_data <- log2(model_data + 1)
  }

  model_data <- model_data %>%
    mutate(across(everything(), as.numeric)) %>%
    drop_na(all_of(target_gene))

  return(model_data)
}

create_learner <- function(model_type, parameters, learner_id = NULL) {
  if (model_type == "glmnet") {
    learner <- lrn("regr.glmnet", alpha = parameters$alpha, s = parameters$s)
  } else if (model_type == "ranger") {
    learner <- lrn("regr.ranger", mtry = parameters$mtry, num.trees = parameters$num.trees, min.node.size = parameters$min.node.size)
  } else if (model_type == "xgboost") {
    learner <- lrn("regr.xgboost", eta = parameters$eta, max_depth = parameters$max_depth, subsample = parameters$subsample, colsample_bytree = parameters$colsample_bytree, lambda = parameters$lambda)
  } else if (model_type == "bart") {
    learner <- lrn("regr.bart", ntree = parameters$ntree, k = parameters$k, power = parameters$power, base = parameters$base)
  } else {
    stop("Unsupported model type: ", model_type)
  }

  if (!is.null(learner_id)) {
    learner$id <- learner_id
  }

  return(learner)
}

run_grid_cv <- function(task, model_type, candidate_grid, folds = 10) {
  resampling <- rsmp("cv", folds = folds)
  resampling$instantiate(task)

  cv_results <- map_dfr(seq_len(nrow(candidate_grid)), function(i) {
    parameters <- as.list(candidate_grid[i, , drop = FALSE])
    learner_id <- paste(model_type, paste(names(parameters), parameters, sep = "_", collapse = "_"), sep = "_")
    learner <- create_learner(model_type, parameters, learner_id)

    resample_result <- resample(task, learner, resampling, store_models = FALSE)

    fold_scores <- as.data.table(resample_result$score(msrs(c("regr.rsq", "regr.mse")))) %>%
      as_tibble() %>%
      transmute(
        model_type = model_type,
        learner_id = learner_id,
        iteration = iteration,
        regr.rsq = regr.rsq,
        regr.mse = regr.mse
      )

    fold_scores <- bind_cols(fold_scores, as_tibble(candidate_grid[rep(i, nrow(fold_scores)), , drop = FALSE]))
    return(fold_scores)
  })

  cv_summary <- cv_results %>%
    group_by(model_type, learner_id) %>%
    summarise(
      mean_rsq = mean(regr.rsq, na.rm = TRUE),
      sd_rsq = sd(regr.rsq, na.rm = TRUE),
      mean_mse = mean(regr.mse, na.rm = TRUE),
      sd_mse = sd(regr.mse, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(candidate_grid %>% mutate(learner_id = map_chr(seq_len(nrow(candidate_grid)), function(i) {
      parameters <- as.list(candidate_grid[i, , drop = FALSE])
      paste(model_type, paste(names(parameters), parameters, sep = "_", collapse = "_"), sep = "_")
    })), by = "learner_id") %>%
    arrange(desc(mean_rsq), mean_mse)

  list(fold_results = cv_results, summary = cv_summary)
}

select_best_learner <- function(task, model_type, cv_summary) {
  best_parameters <- cv_summary %>% slice(1)
  parameter_names <- setdiff(colnames(best_parameters), c("model_type", "learner_id", "mean_rsq", "sd_rsq", "mean_mse", "sd_mse"))
  parameters <- as.list(best_parameters[, parameter_names, drop = FALSE])
  learner <- create_learner(model_type, parameters, learner_id = paste0(model_type, "_best"))
  learner$train(task)
  return(learner)
}

evaluate_prediction <- function(prediction, dataset_name, learner_name) {
  tibble(
    learner = learner_name,
    dataset = dataset_name,
    pearson_r = cor(prediction$response, prediction$truth, method = "pearson", use = "complete.obs"),
    r_squared = prediction$score(msr("regr.rsq")),
    mse = prediction$score(msr("regr.mse"))
  )
}

plot_prediction <- function(prediction, dataset_name, output_file) {
  plot_data <- tibble(
    observed_ERBB2 = prediction$truth,
    predicted_ERBB2 = prediction$response
  )

  prediction_plot <- ggplot(plot_data, aes(x = observed_ERBB2, y = predicted_ERBB2)) +
    geom_point(size = 1.6, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE) +
    stat_cor(method = "pearson") +
    labs(
      title = dataset_name,
      x = "Observed ERBB2 expression",
      y = "Predicted ERBB2 expression"
    ) +
    theme_classic(base_size = 12)

  ggsave(output_file, prediction_plot, width = 4.2, height = 4.2)
  return(prediction_plot)
}

feature_table <- read_csv("input/feature_genes.csv", show_col_types = FALSE)
feature_genes <- unique(feature_table$gene_symbol)

training_expression <- read_expression_matrix("input/training_expression.csv")
external_expression_1 <- read_expression_matrix("input/external_test_expression_1.csv")
external_expression_2 <- read_expression_matrix("input/external_test_expression_2.csv")
external_expression_3 <- read_expression_matrix("input/external_test_expression_3.csv")

training_data <- prepare_model_matrix(training_expression, feature_genes, target_gene)
external_data_1 <- prepare_model_matrix(external_expression_1, feature_genes, target_gene)
external_data_2 <- prepare_model_matrix(external_expression_2, feature_genes, target_gene)
external_data_3 <- prepare_model_matrix(external_expression_3, feature_genes, target_gene)

common_features <- Reduce(intersect, list(
  colnames(training_data),
  colnames(external_data_1),
  colnames(external_data_2),
  colnames(external_data_3)
))

training_data <- training_data[, common_features, drop = FALSE]
external_data_1 <- external_data_1[, common_features, drop = FALSE]
external_data_2 <- external_data_2[, common_features, drop = FALSE]
external_data_3 <- external_data_3[, common_features, drop = FALSE]

training_task <- as_task_regr(training_data, target = target_gene, id = "training_set")
external_task_1 <- as_task_regr(external_data_1, target = target_gene, id = "external_test_set_1")
external_task_2 <- as_task_regr(external_data_2, target = target_gene, id = "external_test_set_2")
external_task_3 <- as_task_regr(external_data_3, target = target_gene, id = "external_test_set_3")

hyperparameter_grids <- list(
  glmnet = expand_grid(
    alpha = c(0, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1),
    s = c(0.0001, 0.001, 0.01, 0.1, 1)
  ),
  ranger = expand_grid(
    mtry = c(3, 6, 10, 20, 40, 60, 80, 95),
    num.trees = c(50, 100, 200, 500, 700, 1000),
    min.node.size = c(1, 3, 5, 7, 10, 15, 20)
  ),
  xgboost = expand_grid(
    eta = c(0.01, 0.1, 0.3, 0.5, 0.7, 1),
    max_depth = c(3, 6, 10, 15, 20),
    subsample = c(0.5, 0.75, 1),
    colsample_bytree = c(0.5, 0.75, 1),
    lambda = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000)
  ),
  bart = expand_grid(
    ntree = c(50, 100, 200, 300, 400, 500),
    k = c(1, 2, 3, 4, 5),
    power = c(0.5, 1, 1.5, 2, 2.5, 3),
    base = c(0.8, 0.85, 0.9, 0.95, 0.99)
  ),
  svm = expand_grid(
    cost = c(1, 5, 10, 50, 100),
    gamma = c(0.0001, 0.001, 0.01, 0.1),
    epsilon = c(0.01, 0.05, 0.1, 0.2)
  )
)

cv_outputs <- imap(hyperparameter_grids, function(candidate_grid, model_type) {
  run_grid_cv(training_task, model_type, candidate_grid, folds = 10)
})

cross_validation_fold_results <- map_dfr(cv_outputs, "fold_results")
cross_validation_summary <- map_dfr(cv_outputs, "summary")

write_csv(cross_validation_fold_results, file.path(output_dir, "cross_validation_fold_results.csv"))
write_csv(cross_validation_summary, file.path(output_dir, "cross_validation_hyperparameter_summary.csv"))

best_learners <- imap(cv_outputs, function(cv_output, model_type) {
  select_best_learner(training_task, model_type, cv_output$summary)
})

ensemble_learner <- gunion(lapply(best_learners, function(learner) po("learner", learner))) %>>%
  po("regravg")
ensemble_learner <- as_learner(ensemble_learner)
ensemble_learner$id <- "mean_ensemble"
ensemble_learner$train(training_task)

all_learners <- c(best_learners, list(mean_ensemble = ensemble_learner))
all_tasks <- list(
  training_set = training_task,
  external_test_set_1 = external_task_1,
  external_test_set_2 = external_task_2,
  external_test_set_3 = external_task_3
)

final_results <- imap_dfr(all_learners, function(learner, learner_name) {
  saveRDS(learner, file.path(model_dir, paste0(learner$id, ".rds")))

  map_dfr(names(all_tasks), function(dataset_name) {
    prediction <- learner$predict(all_tasks[[dataset_name]])

    prediction_table <- tibble(
      sample_id = rownames(all_tasks[[dataset_name]]$data()),
      observed_ERBB2 = prediction$truth,
      predicted_ERBB2 = prediction$response
    )

    write_csv(
      prediction_table,
      file.path(output_dir, paste0(learner$id, "_", dataset_name, "_predictions.csv"))
    )

    plot_prediction(
      prediction,
      paste(learner$id, dataset_name, sep = " - "),
      file.path(output_dir, paste0(learner$id, "_", dataset_name, "_prediction_plot.pdf"))
    )

    evaluate_prediction(prediction, dataset_name, learner$id)
  })
})

write_csv(final_results, file.path(output_dir, "external_validation_summary.csv"))
print(cross_validation_summary %>% group_by(model_type) %>% slice(1) %>% ungroup())
print(final_results)
