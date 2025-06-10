load_all()

my_env <- new.env()
load(paste0(get_methylCIPHER_path(), "/", "SystemsAge_data.RData"), envir = my_env)

debug(calcSystemsAge)

example_SystemsAge <- calcSystemsAge(
  exampleBetas,
  # examplePheno,
  # If `SystemsAge_data.RData` has been downloaded to default folder at $HOME/methylCHIPER
  RData = my_env
)

example_SystemsAge |> class()
example_SystemsAge

pc_env <- new.env()
load(paste0(get_methylCIPHER_path(), "/", "CalcAllPCClocks.RData"), envir = pc_env)

examplePheno$Female <- 1
examplePheno$Age <- examplePheno$age

example_PCClocks <- calcPCClocks(
  exampleBetas,
  examplePheno,
  # If `CalcAllPCClocks.RData` has been downloaded to default folder at $HOME/methylCHIPER
  RData = pc_env
)

example_PCClocks

# Bench mark mean impute
with_missing <- as.matrix(exampleBetas)
set.seed(2025)
ampute <- sample(seq_along(with_missing), size = 10000)
with_missing[ampute] <- NA

meanimpute <- function(x){
  apply(x,2,function(z)ifelse(is.na(z),mean(z,na.rm=T),z))
}

meanimpute1 <- function(x){
  na_indices <- which(is.na(x), arr.ind = TRUE)
  column_means <- colMeans(x, na.rm = TRUE)
  x[na_indices] <- column_means[na_indices[, 2]]
  return(x)
}

all.equal(meanimpute(with_missing), meanimpute1(with_missing))

microbenchmark::microbenchmark(
  meanimpute(with_missing),
  meanimpute1(with_missing)
)
