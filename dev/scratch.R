# Todo
## Imputation CpGs default
## Investigate hangs after error caused by something keeps running even when failed.
## Discuss data vs extdata
## Many data objects have errors and are confusing/redundant: CpGs_GrimAge1 GrimAge1_CpGs.
## Too many CpGs objects/coeff objects. Standardize model params across objects
## Too many R files
## zero CpGs vs imputed values. Add args for imputation
## verbose

# Changes
## Change DNAm to matrix
## Enforce Sample_ID column in pheno. Uniqueness not enforced.
## When pheno is included, automatic alignment are performed
## RData is either a list or a path.
## Please change exampleBetas to matrix. and add Sample_ID, Female, Age for examplePheno

library(FlowSorted.Blood.EPIC)
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

grim1 <- calcGrimAgeV1(betas, phenos)
grim2 <- calcGrimAgeV1_1(betas, phenos)

for(i in setdiff(CalcGrimAgeV1$components, c("Age", "Female"))) {
  print(cor(grim1[[i]], grim2[[i]]))
}

grim1

