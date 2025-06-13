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

orig_clock <- calcPCClocks(test_betas, test_pheno, PCClocks_data)
mod_clock <- calcPCClocks1(test_betas, test_pheno, PCClocks_data)

all.equal(orig_clock, mod_clock)
orig_clock
mod_clock

# git restore --source=a48d583285d92c4f4a1266357dd24a1482894afd -- R/calcPCClocks.R
