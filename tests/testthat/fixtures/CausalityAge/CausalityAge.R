library(tidyverse)
load_all()
GPL21145_pheno <- read_csv(test_path("fixtures", "GPL21145_pheno.csv"))
GPL21145_matrix <- qs2::qs_read(test_path("fixtures", "GPL21145_matrix.qs2"))
scripts <- test_path(
  "fixtures",
  "CausalityAge",
  c("calcCausAge.R", "calcAdaptAge.R", "calcDamAge.R")
)
for (i in scripts) {
  source(i)
}
cpgs <- unique(c(CausAge_CpGs$CpG, AdaptAge_CpGs$CpG, DamAge_CpGs$CpG))
datMeth <- GPL21145_matrix[, intersect(colnames(GPL21145_matrix), cpgs)]
datPheno <- dplyr::select(GPL21145_pheno, Sample_ID = title)

stopifnot(isTRUE(all.equal(row.names(GPL21145_matrix), datPheno$Sample_ID)))
imputeMissingProbes <- vector("numeric", length = length(cpgs))
## `CausAge`
CausAge <- calcCausAge(datMeth, datPheno)
CausAge[, 2] <- CausAge$CausAge
## `AdaptAge`
AdaptAge <- calcAdaptAge(datMeth, datPheno)
AdaptAge[, 2] <- AdaptAge$AdaptAge
## `DamAge`
DamAge <- calcDamAge(datMeth, datPheno)
DamAge[, 2] <- DamAge$DamAge
saveRDS(as.data.frame(CausAge), test_path("fixtures", "CausalityAge", "GPL21145_CausAge.rds"))
saveRDS(as.data.frame(AdaptAge), test_path("fixtures", "CausalityAge", "GPL21145_AdaptAge.rds"))
saveRDS(as.data.frame(DamAge), test_path("fixtures", "CausalityAge", "GPL21145_DamAge.rds"))
qs2::qs_save(datMeth, test_path("fixtures", "CausalityAge", "GPL21145_matrix_CausalityAge.qs2"))
