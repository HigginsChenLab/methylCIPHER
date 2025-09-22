# # TODO
# ### Data and Object Standardization
# - Consolidate data objects for each clock into a single, standardized list with consistent naming (e.g., resolve `CpGs_GrimAge1` vs. `GrimAge1_CpGs`, `CalcGrimAge2` vs. `calcGrimAgeV2`).
# - Standardize model parameters across CpGs and coefficient objects to reduce redundancy.
# - Decide on loading `data` vs. `extdata` for each clock and enforce consistency.
# - Prune the list of exported objects to minimize clutter.
#
# ### Imputation and CpGs Handling
# - Add arguments for CpGs imputation (e.g., zero CpGs vs. imputed values).
# - Set default behavior for CpGs imputation.
# - Add a verbose option to log details (e.g., "Calculating...", "Imputing X CpGs...").
#
# ### Codebase Optimization
# - Reduce the number of R files for better maintainability.
# - Investigate and fix hangs after errors (caused by processes continuing to run post-failure).
#
# ### Discussions
# - Discuss `data` vs. `extdata` usage with the team to finalize approach.
# - duckdb/sql lite to store external data as well as the location to download and imputation standards
# - Setup test suite. Problem: clocks might be inaccurate or dependencies causing things to change
#
# # Changes
# ### Data Input and Alignment
# - Convert `DNAm` to a matrix.
# - Enforce `Sample_ID` column in `pheno` (uniqueness not required).
# - Perform automatic alignment when `pheno` is included.
# - Update `exampleBetas` to a matrix and add `Sample_ID`, `Female`, and `Age` to `examplePheno`.
# - Allow `RData` to be either a list or a file path.
# - Deleted some R files and move functions into one file (trafo and anti.trafo). utils functions into utils.
#
# ### Bug Fixes
# - Fix small bug in `GrimAge1`, and `GrimAge2` where the `temp` object uses `cbind` with `as.numeric`. Specifically, `Female=datPheno$Female` is included but not used in calculations (results unaffected).
#   - Current code: `temp <- as.matrix(cbind(datMeth, Age = as.numeric(datPheno$Age), Female=datPheno$Female))`
#   - Note: `Female` is not used in `CalcGrimAgeV1` or `CalcGrimAge2` model components, so results are correct.
# - Fixed a bug with the pheno[, character] used in subsetting that may cause it to fail with data.table and tibble
# ### Dependency Cleanup
# - Remove `dplyr` and `readr` as dependencies.

load_all()

pc <- load_PCClocks_data()
systems <- load_SystemsAge_data()

calcSystemsAge(DNAm = exampleBetas, pheno = examplePheno, RData = systems)
pheno2 <- examplePheno
pheno2[c(1, 3), "Female"] <- NA
pheno2

calcPCClocks(DNAm = exampleBetas, pheno = pheno2, RData = pc)

calcSystemsAge(DNAm = exampleBetas, ID = "KEKWK", pheno = pheno2, RData = pc)

systems$systems_PCA
