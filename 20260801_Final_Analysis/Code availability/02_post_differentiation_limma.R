# Post-differentiation passage contrasts in log2-normalized expression.
# D0 is excluded; thresholds apply to the equal-weight average contrast.
# Advanced passage denotes the equal-weight average of Middle and Late Passage.

options(stringsAsFactors = FALSE, scipen = 999)

if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
  Sys.setlocale("LC_CTYPE", "Chinese (Traditional)_Taiwan.utf8")
}

required_packages <- c("data.table", "limma")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

# Load Code 1 objects when needed.
code1_objects <- c(
  "input_data",
  "expression_matrix",
  "analysis_expression",
  "sample_metadata",
  "analysis_root"
)

if (!all(vapply(code1_objects, exists, logical(1), inherits = TRUE))) {
  code1_candidates <- c(
    "01_load_data_and_PCA.R",
    file.path("RStudio", "01_load_data_and_PCA.R")
  )
  code1_file <- code1_candidates[file.exists(code1_candidates)][1]
  if (is.na(code1_file)) {
    stop("Run Code 1 first or open this script from the RStudio project")
  }
  source(code1_file)
}

table_directory <- file.path(analysis_root, "Tables")
dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)

# Align metadata and expression columns.
model_metadata <- copy(sample_metadata)
model_metadata <- model_metadata[match(colnames(analysis_expression), sample)]
if (!identical(model_metadata$sample, colnames(analysis_expression))) {
  stop("Sample metadata and expression-matrix columns are not aligned")
}

condition_levels <- as.vector(outer(
  c("Early", "Middle", "Late"),
  c("D0", "D1", "D2", "D6"),
  paste,
  sep = "_"
))
model_metadata[, condition := factor(
  paste(sub(" Passage$", "", as.character(passage)), day, sep = "_"),
  levels = condition_levels
)]

design_matrix <- model.matrix(~0 + condition, data = model_metadata)
colnames(design_matrix) <- levels(model_metadata$condition)

contrast_matrix <- makeContrasts(
  Middle_vs_Early_D1 = Middle_D1 - Early_D1,
  Late_vs_Early_D1 = Late_D1 - Early_D1,
  Middle_vs_Early_D2 = Middle_D2 - Early_D2,
  Late_vs_Early_D2 = Late_D2 - Early_D2,
  Middle_vs_Early_D6 = Middle_D6 - Early_D6,
  Late_vs_Early_D6 = Late_D6 - Early_D6,
  levels = design_matrix
)

base_fit <- lmFit(analysis_expression, design_matrix)
six_contrast_fit <- eBayes(
  contrasts.fit(base_fit, contrast_matrix),
  trend = TRUE,
  robust = TRUE
)

# Joint test of the six post-differentiation contrasts.
joint_results <- as.data.table(
  topTable(
    six_contrast_fit,
    coef = seq_len(ncol(contrast_matrix)),
    number = Inf,
    sort.by = "none"
  ),
  keep.rownames = "ensembl_gene_id"
)
setnames(
  joint_results,
  c("AveExpr", "F", "P.Value", "adj.P.Val"),
  c("average_log2_expression", "joint_F", "joint_p_value", "joint_fdr")
)
joint_results[, (colnames(contrast_matrix)) := NULL]

# Equal-weight average of the six post-differentiation contrasts.
average_contrast <- matrix(
  rowMeans(contrast_matrix),
  ncol = 1,
  dimnames = list(rownames(contrast_matrix), "average_post_effect")
)
average_fit <- eBayes(
  contrasts.fit(base_fit, average_contrast),
  trend = TRUE,
  robust = TRUE
)
average_results <- as.data.table(
  topTable(average_fit, coef = 1, number = Inf, sort.by = "none"),
  keep.rownames = "ensembl_gene_id"
)
average_results <- average_results[, .(
  ensembl_gene_id,
  average_post_log2FC = logFC,
  average_post_t = t,
  average_post_p_value = P.Value,
  average_post_fdr = adj.P.Val
)]

individual_effects <- as.data.table(
  six_contrast_fit$coefficients,
  keep.rownames = "ensembl_gene_id"
)

tested_results <- Reduce(
  function(x, y) merge(x, y, by = "ensembl_gene_id", all = TRUE, sort = FALSE),
  list(joint_results, average_results, individual_effects)
)

effect_columns <- colnames(contrast_matrix)
middle_effect_columns <- paste0("Middle_vs_Early_D", c(1, 2, 6))
late_effect_columns <- paste0("Late_vs_Early_D", c(1, 2, 6))

tested_results[, mean_Middle_vs_Early_log2FC := rowMeans(.SD), .SDcols = middle_effect_columns]
tested_results[, mean_Late_vs_Early_log2FC := rowMeans(.SD), .SDcols = late_effect_columns]
tested_results[, average_fold_change_ratio := 2^average_post_log2FC]

# Two-fold threshold on the average effect.
tested_results[, average_two_fold_pass :=
  average_fold_change_ratio >= 2 |
    average_fold_change_ratio <= 0.5]

tested_results[, positive_contrasts := rowSums(.SD > 0), .SDcols = effect_columns]
tested_results[, negative_contrasts := rowSums(.SD < 0), .SDcols = effect_columns]
tested_results[, trajectory_direction := fcase(
  positive_contrasts == length(effect_columns), "higher",
  negative_contrasts == length(effect_columns), "lower",
  default = "mixed"
)]
tested_results[, all_six_same_direction := trajectory_direction != "mixed"]

tested_results[, joint_selected := joint_fdr < 0.05]
tested_results[, post_selected :=
  average_post_fdr < 0.05 &
    average_two_fold_pass &
    all_six_same_direction]
tested_results[, absolute_average_post_log2FC := abs(average_post_log2FC)]
tested_results[, post_rank := frank(
  average_post_fdr,
  ties.method = "min",
  na.last = "keep"
)]

annotation <- unique(input_data[, .(ensembl_gene_id, symbol)])
annotation[, symbol := fifelse(
  is.na(symbol) | symbol == "",
  ensembl_gene_id,
  symbol
)]

gene_results <- merge(
  annotation,
  tested_results,
  by = "ensembl_gene_id",
  all.x = TRUE,
  sort = FALSE
)
gene_results[, expression_filter_pass := ensembl_gene_id %in% rownames(analysis_expression)]
setorder(
  gene_results,
  -post_selected,
  average_post_fdr,
  -absolute_average_post_log2FC,
  ensembl_gene_id
)

selected_results <- gene_results[post_selected == TRUE]

analysis_summary <- data.table(
  metric = c(
    "input_genes",
    "expression_filter_genes",
    "post_differentiation_contrasts",
    "joint_FDR_cutoff",
    "joint_selected_genes",
    "average_effect_FDR_cutoff",
    "average_fold_change_cutoff",
    "advanced_passage_definition",
    "direction_requirement",
    "post_selected_genes",
    "downregulated_advanced_vs_early_genes",
    "upregulated_advanced_vs_early_genes"
  ),
  value = as.character(c(
    nrow(expression_matrix),
    nrow(analysis_expression),
    ncol(contrast_matrix),
    0.05,
    sum(gene_results$joint_selected, na.rm = TRUE),
    0.05,
    2,
    "equal-weight average of Middle and Late vs Early Passage",
    "all six contrasts must have the same sign",
    nrow(selected_results),
    sum(selected_results$trajectory_direction == "lower"),
    sum(selected_results$trajectory_direction == "higher")
  ))
)

fwrite(
  gene_results,
  file.path(table_directory, "02_post_differentiation_gene_results.csv"),
  bom = TRUE
)
fwrite(
  selected_results,
  file.path(table_directory, "02_selected_consistent_direction_genes.csv"),
  bom = TRUE
)
fwrite(
  analysis_summary,
  file.path(table_directory, "02_differential_analysis_summary.csv"),
  bom = TRUE
)

cat("Code 2 completed successfully.\n")
cat("Expression-filtered genes:", nrow(analysis_expression), "\n")
cat("Joint F-test selected genes:", sum(gene_results$joint_selected, na.rm = TRUE), "\n")
cat("Final consistent-direction genes:", nrow(selected_results), "\n")
cat("  Downregulated advanced vs early:", sum(selected_results$trajectory_direction == "lower"), "\n")
cat("  Upregulated advanced vs early:", sum(selected_results$trajectory_direction == "higher"), "\n")
