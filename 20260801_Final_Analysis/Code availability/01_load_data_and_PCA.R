# Step 1: load the RNA-seq matrix and reproduce the established PCA.
#
# This script is intentionally limited to data loading, sample-level QC and PCA.
# It does not run differential-expression or enrichment analyses.

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260824)

# R launched from some Windows terminals inherits the invalid locale name
# "C.UTF-8".  Set a real UTF-8 Windows locale before resolving this project's
# Chinese path.  RStudio normally already uses a compatible locale.
if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
  Sys.setlocale("LC_CTYPE", "Chinese (Traditional)_Taiwan.utf8")
}

required_packages <- c("readxl", "data.table", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

find_script_path <- function() {
  file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_argument)) {
    return(sub("^--file=", "", file_argument[1]))
  }

  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
    context_path <- tryCatch(
      rstudioapi::getSourceEditorContext()$path,
      error = function(e) ""
    )
    if (nzchar(context_path)) {
      return(context_path)
    }
  }

  ""
}

script_path <- find_script_path()
working_directory <- getwd()
if (basename(working_directory) %in% c("C2C12_RNA-seq analysis", "RStudio")) {
  analysis_root <- dirname(working_directory)
} else if (
  dir.exists(file.path(working_directory, "Figure")) &&
    dir.exists(file.path(working_directory, "Tables"))
) {
  analysis_root <- working_directory
} else if (nzchar(script_path)) {
  analysis_root <- dirname(dirname(script_path))
} else {
  stop("Cannot locate the 20260801_Final_Analysis directory")
}

setwd(analysis_root)
input_file <- file.path("..", "DESeq2_Chiu.xlsx")
if (!file.exists(input_file)) stop("Input file was not found: ", input_file)
figure_directory <- file.path(analysis_root, "Figure")
table_directory <- file.path(analysis_root, "Tables")
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)

cat("Analysis root:", analysis_root, "\n")
cat("Loading:", input_file, "\n")

input_data <- as.data.table(readxl::read_excel(
  input_file,
  sheet = 1,
  .name_repair = "minimal"
))

input_sample_columns <- names(input_data)[grepl("^P", names(input_data))]
if (length(input_sample_columns) != 36L) {
  stop("Expected 36 RNA-seq samples; found ", length(input_sample_columns))
}
if (!all(c("ensembl_gene_id", "symbol") %in% names(input_data))) {
  stop("Required gene-identifier columns are missing")
}
if (anyDuplicated(input_data$ensembl_gene_id)) {
  stop("Duplicated Ensembl gene identifiers were detected")
}

sample_pattern <- "^(P11|P22|P33)_D(0|1|2|6)_([1-3])$"
sample_matches <- regexec(sample_pattern, input_sample_columns)
sample_parts <- regmatches(input_sample_columns, sample_matches)
if (any(lengths(sample_parts) != 4L)) {
  bad_names <- input_sample_columns[lengths(sample_parts) != 4L]
  stop("Unexpected RNA-seq sample names: ", paste(bad_names, collapse = ", "))
}

passage_lookup <- c(
  P11 = "Early Passage",
  P22 = "Middle Passage",
  P33 = "Late Passage"
)
sample_metadata <- data.table(
  passage_code = vapply(sample_parts, `[[`, character(1), 2),
  day_numeric = as.integer(vapply(sample_parts, `[[`, character(1), 3)),
  replicate = as.integer(vapply(sample_parts, `[[`, character(1), 4))
)
sample_metadata[, passage := factor(
  passage_lookup[passage_code],
  levels = c("Early Passage", "Middle Passage", "Late Passage")
)]
sample_metadata[, day := factor(
  paste0("D", day_numeric),
  levels = c("D0", "D1", "D2", "D6")
)]
sample_metadata[, sample := sprintf(
  "%s_D%d_R%02d",
  gsub(" ", "_", as.character(passage)),
  day_numeric,
  replicate
)]
sample_metadata[, passage_code := NULL]
setcolorder(sample_metadata, c("sample", "passage", "day_numeric", "day", "replicate"))
setorder(sample_metadata, passage, day_numeric, replicate)

expression_matrix <- as.matrix(input_data[, ..input_sample_columns])
storage.mode(expression_matrix) <- "double"
rownames(expression_matrix) <- input_data$ensembl_gene_id
colnames(expression_matrix) <- sample_metadata$sample

if (any(!is.finite(expression_matrix)) || any(expression_matrix < 0)) {
  stop("The expression matrix contains invalid values")
}
if (all(expression_matrix == round(expression_matrix))) {
  stop("Integer counts detected; this PCA script expects normalized counts")
}

expression_filter <- rowSums(expression_matrix >= 10) >= 3
log2_expression <- log2(expression_matrix + 1)
analysis_expression <- log2_expression[expression_filter, , drop = FALSE]

gene_variance <- apply(analysis_expression, 1, var)
pca_gene_ids <- names(sort(gene_variance, decreasing = TRUE))
pca_gene_ids <- pca_gene_ids[
  is.finite(gene_variance[pca_gene_ids]) & gene_variance[pca_gene_ids] > 0
]
pca_gene_ids <- head(pca_gene_ids, 1000L)

pca_model <- prcomp(
  t(analysis_expression[pca_gene_ids, , drop = FALSE]),
  center = TRUE,
  scale. = FALSE
)
pca_variance_percent <- 100 * pca_model$sdev^2 / sum(pca_model$sdev^2)

pca_scores <- as.data.table(
  pca_model$x[, 1:5, drop = FALSE],
  keep.rownames = "sample"
)
pca_scores <- merge(sample_metadata, pca_scores, by = "sample", sort = FALSE)
pca_scores[, `:=`(
  total_normalized_count = colSums(expression_matrix)[sample],
  detected_genes = colSums(expression_matrix > 0)[sample],
  median_log2_expression = apply(log2_expression, 2, median)[sample]
)]
setorder(pca_scores, passage, day_numeric, replicate)

# Gene-level audit table for Step 1. This preserves the source gene order and
# makes the expression filter and PCA feature selection directly inspectable.
gene_filter_table <- data.table(
  ensembl_gene_id = input_data$ensembl_gene_id,
  symbol = input_data$symbol,
  mean_normalized_count = rowMeans(expression_matrix),
  median_normalized_count = apply(expression_matrix, 1, median),
  maximum_normalized_count = apply(expression_matrix, 1, max),
  samples_with_normalized_count_ge_10 = rowSums(expression_matrix >= 10),
  expression_filter_pass = as.logical(expression_filter)
)
gene_filter_table[
  is.na(symbol) | symbol == "",
  symbol := ensembl_gene_id
]

variance_rank <- seq_along(variance_order <- names(sort(gene_variance, decreasing = TRUE)))
names(variance_rank) <- variance_order
gene_filter_table[, log2_expression_variance := gene_variance[ensembl_gene_id]]
gene_filter_table[, PCA_variance_rank := as.integer(variance_rank[ensembl_gene_id])]
gene_filter_table[, PCA_top1000_used := ensembl_gene_id %in% pca_gene_ids]

retained_gene_table <- gene_filter_table[expression_filter_pass == TRUE]
pca_gene_table <- gene_filter_table[PCA_top1000_used == TRUE]
setorder(pca_gene_table, PCA_variance_rank)

passage_colours <- c(
  "Early Passage" = "#0072B2",
  "Middle Passage" = "#E69F00",
  "Late Passage" = "#CC79A7"
)
day_shapes <- c(D0 = 21, D1 = 22, D2 = 23, D6 = 24)

pca_figure <- ggplot(pca_scores, aes(PC1, PC2, colour = passage, shape = day)) +
  geom_point(size = 3, stroke = 0.8, fill = "white") +
  scale_colour_manual(values = passage_colours) +
  scale_shape_manual(values = day_shapes) +
  labs(
    x = sprintf("PC1 (%.1f%%)", pca_variance_percent[1]),
    y = sprintf("PC2 (%.1f%%)", pca_variance_percent[2]),
    colour = "Passage",
    shape = "Day"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11, face = "bold")
  )

pca_figure_file <- file.path(figure_directory, "01_RNA_PCA.png")
ggsave(
  pca_figure_file,
  plot = pca_figure,
  width = 7.4,
  height = 5.2,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(figure_directory, "01_RNA_PCA.pdf"),
  plot = pca_figure,
  width = 7.4,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)

loading_summary <- data.table(
  metric = c(
    "input_file",
    "input_md5",
    "input_genes",
    "samples",
    "expression_filter_rule",
    "expression_filter_genes",
    "PCA_genes",
    "PCA_centered",
    "PCA_scaled",
    "PC1_variance_percent",
    "PC2_variance_percent"
  ),
  value = c(
    input_file,
    unname(tools::md5sum(input_file)),
    nrow(expression_matrix),
    ncol(expression_matrix),
    "normalized count >= 10 in >= 3 samples",
    nrow(analysis_expression),
    length(pca_gene_ids),
    "TRUE",
    "FALSE",
    sprintf("%.6f", pca_variance_percent[1]),
    sprintf("%.6f", pca_variance_percent[2])
  )
)

fwrite(
  sample_metadata,
  file.path(table_directory, "01_sample_metadata.csv"),
  bom = TRUE
)
fwrite(
  gene_filter_table,
  file.path(table_directory, "01_gene_filter_all_genes.csv"),
  bom = TRUE
)
fwrite(
  retained_gene_table,
  file.path(table_directory, "01_expression_filter_retained_genes.csv"),
  bom = TRUE
)
fwrite(
  pca_gene_table,
  file.path(table_directory, "01_PCA_top1000_gene_table.csv"),
  bom = TRUE
)
fwrite(
  pca_scores,
  file.path(table_directory, "01_sample_QC_and_PCA_scores.csv"),
  bom = TRUE
)
fwrite(
  loading_summary,
  file.path(table_directory, "01_data_loading_and_PCA_summary.csv"),
  bom = TRUE
)

cat("Loaded genes:", nrow(expression_matrix), "\n")
cat("Samples:", ncol(expression_matrix), "\n")
cat("Genes after expression filter:", nrow(analysis_expression), "\n")
cat("PCA genes:", length(pca_gene_ids), "\n")
cat(sprintf("PC1: %.1f%%; PC2: %.1f%%\n", pca_variance_percent[1], pca_variance_percent[2]))
cat("PCA figure:", pca_figure_file, "\n")
cat("Step 1 completed successfully.\n")
