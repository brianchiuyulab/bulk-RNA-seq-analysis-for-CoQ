# C2C12 RNA-seq analysis

Open `RStudio/C2C12_RNA-seq analysis.Rproj` and run `00_run_all.R`, or run Codes 01–03 in order.

## Analysis workflow

1. `01_load_data_and_PCA.R` loads normalized expression, applies the expression filter and performs PCA using the 1,000 genes with the highest variance.
2. `02_post_differentiation_limma.R` fits limma models to log2 normalized expression for Middle-versus-Early and Late-versus-Early contrasts at D1, D2 and D6. Final genes require BH-FDR below 0.05 for the equal-weight average contrast, an absolute average change of at least two-fold and the same direction across all six contrasts.
3. `03_global_heatmap_KEGG.R` draws the 1,498-gene heatmap and performs separate KEGG analyses for the higher and lower signatures using all expression-filtered genes as the background.

Passage groups are reported as Early Passage, Middle Passage and Late Passage throughout the processed data, contrasts and tables. Dense figure labels use `Early P.`, `Middle P.` and `Late P.`. Axis, legend and annotation text is enlarged and bold.

## Directories

- `RStudio`: working scripts and RStudio project.
- `Code availability`: manuscript-ready copies of the same scripts.
- `Figure`: complete PNG and PDF outputs.
- `Tables`: complete CSV outputs.
