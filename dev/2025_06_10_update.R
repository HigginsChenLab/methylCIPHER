# Report ----
remove.packages("methylCIPHER")
remotes::install_github("hhp94/methylCIPHER", ref = "mod")
library(methylCIPHER)

# mean impute vectorization ----
with_missing <- as.matrix(exampleBetas)
set.seed(2025)
ampute <- sample(seq_along(with_missing), size = 1000) # inject more NA value
with_missing[ampute] <- NA

## old
meanimpute <- function(x){
  apply(x,2,function(z)ifelse(is.na(z),mean(z,na.rm=T),z))
}

## new vectorized
meanimpute1 <- function(x){
  na_indices <- which(is.na(x), arr.ind = TRUE)
  column_means <- colMeans(x, na.rm = TRUE)
  x[na_indices] <- column_means[na_indices[, 2]]
  return(x)
}

## same results
all.equal(meanimpute(with_missing), meanimpute1(with_missing))

## Bench mark mean impute. New implementation 10 - 100x faster
microbenchmark::microbenchmark(
  meanimpute(with_missing),
  meanimpute1(with_missing)
)

# manual download weights ----
## Added support for manual pathing. Download the `SystemsAge_data.RData` or
## `CalcAllPCClocks.RData` and add either a path, list, or environment to the
## `RData` argument

# Default path to download to can be calculated here
get_methylCIPHER_default_path()

# Let's say you downloaded the weights to the current folder, you can change as follows
options(methylCIPHER.path = ".")
get_methylCIPHER_path() # This shows the current set path. If not set then is default.

v1 <- calcSystemsAge(
  exampleBetas,
  RData = get_methylCIPHER_default_path() # <----- pass the path where the files are downloaded here
)

# pass an environment instead
my_env <- new.env()
load(
  paste0(
    get_methylCIPHER_default_path(),
    "/",
    "SystemsAge_data.RData"
    ),
  envir = my_env
)
names(my_env)

v2 <- calcSystemsAge(
  exampleBetas,
  RData = my_env # <----- passing in the environment prevents loading the object again
)

# pass weight as a list also works. This is important because we will not use environment
# in the future
load(
  paste0(
    get_methylCIPHER_default_path(),
    "/",
    "SystemsAge_data.RData"
    )
)

weight_list <- list(
  CpGs = CpGs,
  meanimpute = meanimpute,
  Systems_clock_coefficients = Systems_clock_coefficients,
  DNAmPCA = DNAmPCA,
  Age_prediction_model = Age_prediction_model,
  transformation_coefs = transformation_coefs,
  Predicted_age_coefficients = Predicted_age_coefficients,
  imputeMissingCpGs = imputeMissingCpGs,
  system_vector_coefficients = system_vector_coefficients,
  system_scores_coefficients_scale = system_scores_coefficients_scale,
  systems_PCA = systems_PCA
)

v3 <- calcSystemsAge(
  exampleBetas,
  RData = weight_list # <----- passing in the list prevents loading the object again
)

all.equal(v1, v2)
all.equal(v1, v3)
