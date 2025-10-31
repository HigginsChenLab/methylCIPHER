library(data.table)
library(glue)

fixtures <- "tests/testthat/fixtures"
# The matrix will fail to download. This is just for the pheno data
if (!file.exists(glue("{fixtures}/GPL21145.soft.gz"))) {
  library(GEOquery)

  gse <- getGEO(
    "GSE286313",
    destdir = glue("{fixtures}"),
    GSEMatrix = TRUE,
    getGPL = TRUE,
    parseCharacteristics = TRUE
  )

  GPL21145_pheno <- pData(gse$`GSE286313-GPL21145_series_matrix.txt.gz`)
  data.table::fwrite(GPL21145_pheno, glue("{fixtures}/GPL21145_pheno.csv"))
}

# EPICv1
EPICv1 <- glue("{fixtures}/GSE286313_MatrixProcessed_GPL21145.csv")
if (!file.exists(EPICv1)) {
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286313/suppl/GSE286313%5FMatrixProcessed%5FGPL21145%2Ecsv%2Egz",
    destfile = EPICv1
  )
  gunzip(EPICv1)
}

# Keep all CpGs, don't do Pval gating
GPL21145 <- fread(EPICv1)
GPL21145_pheno <- fread(glue("{fixtures}/GPL21145_pheno.csv"))
GPL21145_matrix <- as.matrix(GPL21145[, .SD, .SDcols = GPL21145_pheno$title])
row.names(GPL21145_matrix) <- GPL21145$ID_REF
GPL21145_matrix <- t(GPL21145_matrix)
qs2::qs_save(GPL21145_matrix, glue("{fixtures}/GPL21145_matrix.qs2"))
GPL21145_matrix[1:5, 1:5]
