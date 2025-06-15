# Todo
## Imputation CpGs default
## Investigate hangs after error caused by something keeps running even when failed.
## Discuss data vs extdata
## Many data objects have errors and are confusing/redundant: CpGs_GrimAge1 GrimAge1_CpGs. `CalcGrimAge2` is different than `calcGrimAgeV2`, etc.
### Each clock should have just one object as dependency. This should be a list. This should have standardized names. And load extdata or not
## Too many CpGs objects/coeff objects. Standardize model params across objects
## Too many R files
## zero CpGs vs imputed values. Add args for imputation
## verbose city option? Calculating? Imputing X CpGs?
## Prune list of export

# Changes
## Change DNAm to matrix
## Enforce Sample_ID column in pheno. Uniqueness not enforced.
## When pheno is included, automatic alignment are performed
## RData is either a list or a path.
## Please change exampleBetas to matrix. and add Sample_ID, Female, Age for examplePheno
## Small bug with the "temp" object of PCCClocks, SystemsAge, GrimAge1, and GrimAge2 of cbind-ing with as.numeric
## Female is not included in the temp object, but since the clocks don't use Female in the calculation of the
## clocks, it doesn't affect theresults
## remove dplyr and readr as dependencies

## In here, Female=datPheno$Female is dropped from as.numeric.
## temp <- as.matrix(cbind(datMeth, Age = as.numeric(datPheno$Age, Female=datPheno$Female)))

## But since we see Female is in non of the components, the results are still correct
# lapply(
#   CalcGrimAgeV1[grepv("model", names(CalcGrimAgeV1))],
#   \(x) { names(x)[grepl("Age|Female", names(x))]}
# )

## But since we see Female is in non of the components, the results are still correct
# lapply(
#   CalcGrimAge2[grepv("model", names(CalcGrimAge2))],
#   \(x) { names(x)[grepl("Age|Female", names(x))]}
# )

library(FlowSorted.Blood.EPIC)
library(data.table)

FlowSorted.Blood.EPIC <- libraryDataGet("FlowSorted.Blood.EPIC")

load_all()

phenos <- as.data.frame(colData(FlowSorted.Blood.EPIC))
phenos <- subset(phenos, select = c("Sex", "Age"))
phenos$Sample_ID <- row.names(phenos)
phenos$Female <- as.numeric(phenos$Sex == "F")
row.names(phenos) <- NULL
phenos <- subset(phenos, !is.na(Sex), select = c("Age", "Sample_ID", "Female"))

betas <- getBeta(FlowSorted.Blood.EPIC)
betas <- removeNAcol(t(betas)[phenos$Sample_ID, , drop = F])

betas |> dim()

library(tibble)
grim1 <- calcGrimAgeV2(betas, phenos)
grim2 <- calcGrimAgeV2_1(betas, as.data.table(phenos))

for(i in setdiff(CalcGrimAge2$components, c("Age", "Female"))) {
  print(cor(grim1[[i]], grim2[[i]]))
}

plot(grim1[[i]], grim2[[i]])
abline(a = 0, b = 1)

# # formatHorvath
# formatHorvathOnline(exampleBetas, examplePheno, ".")
# formatHorvathOnline_1(exampleBetas, examplePheno, "new")
#
# rlang::hash_file("DNAm_Horvath_Online_Calculator_input.csv") == rlang::hash_file("new/DNAm_Horvath_Online_Calculator_input.csv")
# rlang::hash_file("Pheno_Horvath_Online_Calculator_input.csv") == rlang::hash_file("new/Pheno_Horvath_Online_Calculator_input.csv")
