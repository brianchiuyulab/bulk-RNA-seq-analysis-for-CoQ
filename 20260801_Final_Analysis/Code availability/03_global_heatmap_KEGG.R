# Global heatmap and KEGG over-representation analysis.

options(stringsAsFactors = FALSE, scipen = 999)

if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
  Sys.setlocale("LC_CTYPE", "Chinese (Traditional)_Taiwan.utf8")
}

required_packages <- c(
  "data.table", "ComplexHeatmap", "circlize", "clusterProfiler",
  "AnnotationDbi", "org.Mm.eg.db", "ggplot2"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(ggplot2)
})

# Load Code 2 objects when needed.
code2_objects <- c(
  "input_data", "analysis_expression", "model_metadata",
  "gene_results", "selected_results", "analysis_root"
)
if (!all(vapply(code2_objects, exists, logical(1), inherits = TRUE))) {
  code2_candidates <- c(
    "02_post_differentiation_limma.R",
    file.path("RStudio", "02_post_differentiation_limma.R")
  )
  code2_file <- code2_candidates[file.exists(code2_candidates)][1]
  if (is.na(code2_file)) stop("Run Code 2 first or open this script from the RStudio project")
  source(code2_file)
}

figure_directory <- file.path(analysis_root, "Figure")
table_directory <- file.path(analysis_root, "Tables")
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)

# Heatmap of all genes selected in Code 2.
selected_gene_ids <- intersect(
  selected_results$ensembl_gene_id,
  rownames(analysis_expression)
)
heatmap_matrix <- analysis_expression[selected_gene_ids, , drop = FALSE]
heatmap_matrix <- t(scale(t(heatmap_matrix)))
heatmap_matrix <- heatmap_matrix[complete.cases(heatmap_matrix), , drop = FALSE]
heatmap_matrix[heatmap_matrix > 2] <- 2
heatmap_matrix[heatmap_matrix < -2] <- -2

heatmap_gene_table <- selected_results[
  match(rownames(heatmap_matrix), ensembl_gene_id),
  .(
    ensembl_gene_id, symbol, trajectory_direction,
    average_post_log2FC, average_post_fdr
  )
]
fwrite(
  heatmap_gene_table,
  file.path(table_directory, "03_heatmap_gene_table.csv"),
  bom = TRUE
)

heatmap_metadata <- copy(model_metadata)
setorder(heatmap_metadata, passage, day_numeric, replicate)
heatmap_matrix <- heatmap_matrix[, heatmap_metadata$sample, drop = FALSE]
heatmap_metadata[, passage_short := factor(
  sub(" Passage$", " P.", as.character(passage)),
  levels = c("Early P.", "Middle P.", "Late P.")
)]
heatmap_metadata[, figure_sample := sprintf(
  "%s · %s · R%02d", passage_short, day, replicate
)]

passage_colors <- c(
  "Early Passage" = "#3B6FB6",
  "Middle Passage" = "#E69F00",
  "Late Passage" = "#B24C63"
)
day_colors <- c(D0 = "#D9D9D9", D1 = "#7FC97F", D2 = "#BEAED4", D6 = "#FDCDAC")
column_annotation <- HeatmapAnnotation(
  Passage = heatmap_metadata$passage,
  Day = heatmap_metadata$day,
  col = list(Passage = passage_colors, Day = day_colors),
  annotation_name_gp = grid::gpar(fontsize = 11, fontface = "bold"),
  annotation_legend_param = list(
    Passage = list(title_gp = grid::gpar(fontsize = 11, fontface = "bold"), labels_gp = grid::gpar(fontsize = 10, fontface = "bold")),
    Day = list(title_gp = grid::gpar(fontsize = 11, fontface = "bold"), labels_gp = grid::gpar(fontsize = 10, fontface = "bold"))
  )
)

heatmap_object <- Heatmap(
  heatmap_matrix,
  name = "Row Z-score",
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B")),
  top_annotation = column_annotation,
  column_split = heatmap_metadata$passage,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_labels = heatmap_metadata$figure_sample,
  column_names_gp = grid::gpar(fontsize = 8, fontface = "bold"),
  column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
  heatmap_legend_param = list(
    title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 10, fontface = "bold")
  ),
  use_raster = TRUE,
  raster_quality = 3,
  border = FALSE
)

png(
  file.path(figure_directory, "03_all_selected_genes_heatmap.png"),
  width = 14, height = 10, units = "in", res = 300, type = "cairo-png"
)
draw(heatmap_object, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

pdf(
  file.path(figure_directory, "03_all_selected_genes_heatmap.pdf"),
  width = 14, height = 10, useDingbats = FALSE
)
draw(heatmap_object, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# KEGG analysis with all tested genes as the background.
tested_gene_ids <- rownames(analysis_expression)
id_map <- as.data.table(AnnotationDbi::select(
  org.Mm.eg.db,
  keys = tested_gene_ids,
  keytype = "ENSEMBL",
  columns = "ENTREZID"
))
setnames(id_map, c("ENSEMBL", "ENTREZID"), c("ensembl_gene_id", "entrez_id"))
id_map <- unique(id_map[!is.na(entrez_id) & entrez_id != ""])

selected_map <- merge(
  selected_results[, .(ensembl_gene_id, trajectory_direction)],
  id_map,
  by = "ensembl_gene_id",
  allow.cartesian = TRUE
)
background_entrez <- unique(id_map$entrez_id)
higher_label <- "Upregulated advanced-passage signature"
lower_label <- "Downregulated advanced-passage signature"
gene_sets <- setNames(
  list(
    unique(selected_map[trajectory_direction == "higher", entrez_id]),
    unique(selected_map[trajectory_direction == "lower", entrez_id])
  ),
  c(higher_label, lower_label)
)

run_kegg <- function(entrez_ids, gene_set) {
  result <- enrichKEGG(
    gene = entrez_ids,
    organism = "mmu",
    keyType = "ncbi-geneid",
    universe = background_entrez,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500
  )
  result <- as.data.table(as.data.frame(result))
  if (nrow(result)) result[, gene_set := gene_set]
  result
}

kegg_results <- rbindlist(
  Map(run_kegg, gene_sets, names(gene_sets)),
  fill = TRUE,
  use.names = TRUE
)
setorder(kegg_results, gene_set, p.adjust, pvalue, -Count)

fwrite(
  kegg_results[gene_set == higher_label],
  file.path(table_directory, "03_KEGG_upregulated_advanced_vs_early.csv"),
  bom = TRUE
)
fwrite(
  kegg_results[gene_set == lower_label],
  file.path(table_directory, "03_KEGG_downregulated_advanced_vs_early.csv"),
  bom = TRUE
)

top_kegg <- kegg_results[p.adjust < 0.05, head(.SD, 15), by = gene_set]
fwrite(top_kegg, file.path(table_directory, "03_KEGG_top_pathways.csv"), bom = TRUE)
if (nrow(top_kegg)) {
  top_kegg[, gene_ratio_numeric := vapply(
    strsplit(GeneRatio, "/", fixed = TRUE),
    function(x) as.numeric(x[1]) / as.numeric(x[2]),
    numeric(1)
  )]
  top_kegg[, pathway_display := sub(
    " - Mus musculus \\(house mouse\\)$", "", Description
  )]
  top_kegg[, pathway_label := paste(pathway_display, gene_set, sep = "___")]
  top_kegg[, pathway_label := factor(pathway_label, levels = rev(unique(pathway_label)))]

  kegg_figure <- ggplot(
    top_kegg,
    aes(x = gene_ratio_numeric, y = pathway_label, size = Count, color = -log10(p.adjust))
  ) +
    geom_point(alpha = 0.9) +
    facet_grid(gene_set ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_color_viridis_c(option = "C", name = expression(-log[10](FDR))) +
    scale_size_continuous(name = "Gene count") +
    labs(x = "Gene ratio", y = NULL) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.major.y = element_blank(),
      strip.background = element_rect(fill = "grey92", color = "grey70"),
      strip.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11, face = "bold")
    )

  ggsave(
    file.path(figure_directory, "03_KEGG_top_pathways.png"),
    kegg_figure, width = 10, height = 8, units = "in", dpi = 300, bg = "white"
  )
  ggsave(
    file.path(figure_directory, "03_KEGG_top_pathways.pdf"),
    kegg_figure, width = 10, height = 8, units = "in", device = cairo_pdf
  )
}

analysis_summary <- data.table(
  metric = c(
    "heatmap_genes", "heatmap_samples", "tested_genes",
    "background_entrez_ids", "selected_entrez_ids",
    "advanced_passage_definition",
    "upregulated_advanced_vs_early_entrez_ids",
    "downregulated_advanced_vs_early_entrez_ids",
    "KEGG_FDR_cutoff",
    "significant_KEGG_upregulated_advanced_vs_early",
    "significant_KEGG_downregulated_advanced_vs_early"
  ),
  value = as.character(c(
    nrow(heatmap_matrix), ncol(heatmap_matrix), length(tested_gene_ids),
    length(background_entrez), length(unique(selected_map$entrez_id)),
    "equal-weight average of Middle and Late vs Early Passage",
    length(gene_sets[[higher_label]]), length(gene_sets[[lower_label]]),
    0.05,
    nrow(kegg_results[gene_set == higher_label & p.adjust < 0.05]),
    nrow(kegg_results[gene_set == lower_label & p.adjust < 0.05])
  ))
)
fwrite(analysis_summary, file.path(table_directory, "03_global_analysis_summary.csv"), bom = TRUE)

cat("Code 3 completed successfully.\n")
cat("Heatmap genes:", nrow(heatmap_matrix), "\n")
cat("Significant KEGG pathways (upregulated/downregulated advanced vs early):",
    nrow(kegg_results[gene_set == higher_label & p.adjust < 0.05]), "/",
    nrow(kegg_results[gene_set == lower_label & p.adjust < 0.05]), "\n")
