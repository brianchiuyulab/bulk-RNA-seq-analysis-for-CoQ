scripts <- c(
  "01_load_data_and_PCA.R",
  "02_post_differentiation_limma.R",
  "03_global_heatmap_KEGG.R"
)

script_directory <- getwd()
if (!all(file.exists(scripts))) stop("Open the RNA-seq RStudio project")
for (script in scripts) {
  setwd(script_directory)
  source(script, local = .GlobalEnv)
}

cat("RNA-seq analysis completed successfully.\n")
