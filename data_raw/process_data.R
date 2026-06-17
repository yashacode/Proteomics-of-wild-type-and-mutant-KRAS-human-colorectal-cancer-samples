setwd('data_raw')

library(dplyr)

prot_data <- data.table::fread("report1.pg_matrix.tsv")[,1:32]

# Extract sample ID (M######a) from full file path column names
colnames(prot_data) <- vapply(colnames(prot_data), function(s) {
  m <- regmatches(s, regexpr("M[0-9]+a", s))
  if (length(m) > 0) m else s
}, character(1))

prot_data[, Protein.Group := sub(";.*", "", Protein.Group)]

metadata <- as.data.frame(data.table::fread("sdrf-human.tsv"))
meta <- data.table::data.table(
  SampleID = vapply(metadata$`source name`, function(s){
    parts <- strsplit(s, "-", fixed = TRUE)[[1]]
    if (length(parts) >= 3) parts[2] else NA_character_
  }, character(1)),

  mutation = vapply(metadata$`source name`, function(s) {
    parts <- strsplit(s, "-", fixed = TRUE)[[1]]
    if (length(parts) >= 3) paste(parts[3:length(parts)], collapse = "-") else NA_character_
  }, character(1))
)
meta$mutated <- ifelse(
  trimws(meta$mutation) == "KRAS/Wildtype(Tumor)",
  "Wildtype",
  "mutant"
)

se <- ProtPipe::create_se(as.data.frame(prot_data), as.data.frame(meta))

# impute_left_dist() draws random values, so seed for reproducible results
set.seed(42)
se_p <- se %>%
  ProtPipe::filter_proteins_by_percent(50) %>%
  ProtPipe::impute_left_dist()

# Save the heatmap inputs as plain objects (matrix + annotations) rather
# than a SummarizedExperiment, so the deployed app needs no Bioconductor
# dependency. Components stay row/column aligned.
heatmap_data <- list(
  expr     = as.matrix(SummarizedExperiment::assay(se_p)),
  row_data = as.data.frame(SummarizedExperiment::rowData(se_p)),
  col_data = as.data.frame(SummarizedExperiment::colData(se_p))
)
saveRDS(heatmap_data, "../inst/app/data/heatmap_data.rds")

DE <- ProtPipe::do_limma_binary(se_p, "mutated", "mutant", "Wildtype")

# Keep the full DE result, including the per-sample intensity columns
# (M#####a). The volcano only uses logFC / P.Value / Genes, but the
# sample columns drive the dataset stat boxes and the Feature Explorer.
cancer_dataset <- DE

saveRDS(cancer_dataset, "../inst/app/data/cancer_de.rds")

# Sample -> group mapping, used by the protein-detail boxplot to split
# per-sample intensities into mutant vs wildtype. SampleID matches the
# M#####a intensity column names in the DE table.
sample_meta <- data.frame(
  SampleID = meta$SampleID,
  group = ifelse(meta$mutated == "Wildtype", "Wildtype", "Mutant"),
  stringsAsFactors = FALSE
)
saveRDS(sample_meta, "../inst/app/data/sample_meta.rds")
