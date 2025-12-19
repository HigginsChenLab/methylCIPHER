library(tidyverse)
load_all()

GPL21145_pheno <- read_csv(test_path("fixtures", "GPL21145_pheno.csv"))
GPL21145_matrix <- qs2::qs_read(test_path("fixtures", "GPL21145_matrix.qs2"))

scripts <- test_path(
  "fixtures",
  "SenescenceAge",
  c("calcSenChronoAge.R", "calcSenCultureAge.R", "calcSenMortalityAge.R")
)

for (i in scripts) {
  source(i)
}

cpgs <- unique(
  c(
    row.names(SenChronoAge_CpGs),
    row.names(SenCultureAge_CpGs),
    row.names(SenMortalityAge_CpGs)
  )
)

datMeth <- GPL21145_matrix[, intersect(colnames(GPL21145_matrix), cpgs)]
datPheno <- select(GPL21145_pheno, Sample_ID = title)
stopifnot(isTRUE(all.equal(row.names(GPL21145_matrix), datPheno$Sample_ID)))

# The original scripts are very bugged. The `imputeMissingProbes` object is missing
# `(Intercept)`. The function is also using the wrong column name.
imputeMissingProbes <- vector("numeric", length = length(cpgs))
names(SenChronoAge_CpGs)[1] <- "coefficient"
names(SenCultureAge_CpGs)[1] <- "coefficient"
names(SenMortalityAge_CpGs)[1] <- "coefficient"

## `SenChronoAge`
SenChronoAge <- calcSenChronoAge(datMeth, datPheno)
SenChronoAge[, 2] <- SenChronoAge$SenChronoAge
## `SenMortalityAge`
SenMortalityAge <- calcSenMortalityAge(datMeth, datPheno)
SenMortalityAge[, 2] <- SenMortalityAge$SenMortalityAge
## `SenCultureAge`
SenCultureAge <- calcSenCultureAge(datMeth, datPheno)
SenCultureAge[, 2] <- SenCultureAge$SenCultureAge

saveRDS(as.data.frame(SenChronoAge), test_path("fixtures", "SenescenceAge", "GPL21145_SenChronoAge.rds"))
saveRDS(as.data.frame(SenMortalityAge), test_path("fixtures", "SenescenceAge", "GPL21145_SenMortalityAge.rds"))
saveRDS(as.data.frame(SenCultureAge), test_path("fixtures", "SenescenceAge", "GPL21145_SenCultureAge.rds"))
qs2::qs_save(datMeth, test_path("fixtures", "SenescenceAge", "GPL21145_matrix_SenescenceAge.qs2"))
