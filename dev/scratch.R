# Todo
## Imputation CpGs default
## Investigate hangs after error caused by something keeps running even when failed.

# Changes
## Change DNAm to matrix
## Enforce Sample_ID column in pheno. Uniqueness not enforced.
## When pheno is included, automatic alignment are performed
## RData is either a list or a path.
## Please change exampleBetas to matrix. and add Sample_ID, Female, Age for examplePheno

load_all()
PCClocks_data <- load_PCClocks_data()

# original ----
test_betas <- impute_DNAm(exampleBetas, method = "mean", CpGs = PCClocks_data$imputeMissingCpGs)
test_pheno <- align_pheno(examplePheno, row.names(test_betas))

mod_clock <- calcPCClocks(test_betas, test_pheno, PCClocks_data)

