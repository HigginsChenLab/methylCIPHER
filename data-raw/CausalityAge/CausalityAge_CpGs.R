library(data.table)

prep <- function(path) {
  dt <- fread(path)
  setnames(dt, c("CpG", "Beta"))
  intercept <- dt[CpG == "(Intercept)", Beta]
  dt <- dt[CpG != "(Intercept)"]
  list(cpgs = as.data.frame(dt), intercept = intercept)
}

caus <- prep("data-raw/CausalityAge/YingCausAge.csv")
adapt <- prep("data-raw/CausalityAge/YingAdaptAge.csv")
dam <- prep("data-raw/CausalityAge/YingDamAge.csv")

CausAge_CpGs <- caus$cpgs
AdaptAge_CpGs <- adapt$cpgs
DamAge_CpGs <- dam$cpgs

message("CausAge intercept:  ", caus$intercept)
message("AdaptAge intercept: ", adapt$intercept)
message("DamAge intercept:   ", dam$intercept)

use_data(CausAge_CpGs, AdaptAge_CpGs, DamAge_CpGs, overwrite = TRUE)
