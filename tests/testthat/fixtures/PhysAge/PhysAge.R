# Create DNAm Physiological health Age (PhysAge) and DNAm surrogate scores
# Em Arpawong - 10/2/25

################################################
# Recommended steps:
# Step 0: Median impute CpGs that are missing for some samples
# Step 1: Impute CpGs missing entirely from your dataset using the "HRS_cpg_training_means.txt" provided
# Step 2: Calculate DNAm raw surrogate scores
# Step 3: Calculate continuous sum of z-scores for DNAmPhysAge score & convert to metric of years

################################################

# Note: There are 2 places below that require your input, noted by asterisks ******

rm(list = ls())
library(dplyr)
library(tibble)
library(purrr)

# code assumes you have all data and weight files in the same directory

# load your DNAm data with rows as samples (samples as rownames) and columns as probe IDs (colnames as CpGs)
betas <- qs2::qs_read(test_path("fixtures", "GPL21145_matrix.qs2"))

# limit to the target CpGs needed for efficiency
want <- read.table("data-raw/PhysAge/CpG lists and weights/HRS_cpg_training_means.txt", header = TRUE)$cpg
matched <- intersect(colnames(betas), want)
x <- betas[, matched]

################################################

weight_files_names <- c(
  DNAmPeakflow   = "peakflow_DNAm_weights.txt",
  DNAmHbA1c      = "hba1c_DNAm_weights.txt",
  DNAmHDL        = "PHDLD_DNAm_weights.txt",
  DNAmPulsePr    = "pulsepressure_DNAm_weights.txt",
  DNAmCRP        = "PCRP_v3_DNAm_weights.txt",
  DNAmCystatinC  = "PCYSC_v2_DNAm_weights.txt",
  DNAmDHEAS      = "PDHEASE_DNAm_weights.txt",
  DNAmWHR        = "WHR_DNAm_weights.txt"
)

weight_files <- paste0("data-raw/PhysAge/CpG lists and weights", "/", weight_files_names)
names(weight_files) <- names(weight_files_names)

scores_all <- imap(weight_files, function(file, label) {
  wf <- read.table(file, col.names = c("cpgs", "weights"))
  intercept <- as.numeric(wf$weights[1])
  w <- wf[-1, c("cpgs", "weights")]

  # keep only CpGs present in x, and align weights to those columns
  keep <- intersect(colnames(x), w$cpgs)
  if (length(keep) == 0L) {
    # no overlapping CpGs → return NAs for this score
    return(tibble(
      sample = rownames(x),
      !!paste0(label, "_final") := NA_real_,
      !!paste0(label, "_raw") := NA_real_
    ))
  }

  w_aligned <- w$weights[match(keep, w$cpgs)]
  # multiply each CpG column by its weight, then take (non-NA) row means
  raw <- rowMeans(sweep(x[, keep, drop = FALSE], 2, w_aligned, `*`), na.rm = TRUE)

  tibble(
    sample = rownames(x),
    !!paste0(label, "_final") := raw + intercept,
    !!paste0(label, "_raw") := raw
  )
}) %>%
  reduce(full_join, by = "sample") %>%
  column_to_rownames("sample") %>%
  .[order(rownames(.)), , drop = FALSE]

# reverse-scale and composite, unchanged (minor tidy tweak)
neg_vars <- c("DNAmPeakflow_raw", "DNAmDHEAS_raw", "DNAmHDL_raw")

scores_all <- scores_all %>%
  mutate(across(
    ends_with("_raw"),
    ~ as.vector(scale(.)) * if_else(cur_column() %in% neg_vars, -1, 1),
    .names = "{paste0('z', sub('_raw$', '', .col))}"
  )) %>%
  mutate(
    DNAmPhysAge  = rowSums(across(starts_with("zDNAm")), na.rm = TRUE),
    zDNAmPhysAge = as.vector(scale(DNAmPhysAge))
  ) %>%
  mutate(across(
    starts_with("zDNAm"),
    ~ (.x * 9.630436) + 68.14726,
    .names = "{sub('^z', '', .col)}_years"
  ))

scores_all$Sample_ID <- row.names(scores_all)
row.names(scores_all) <- NULL

saveRDS(scores_all, test_path("fixtures", "PhysAge", "GPL21145_PhysAge.rds"))
qs2::qs_save(x, test_path("fixtures", "PhysAge", "GPL21145_matrix_PhysAge.qs2"))

# GPL21145_pheno <- data.table::fread(test_path("fixtures", "GPL21145_pheno.csv"))
# GPL21145_PhysAge <- inner_join(scores_all, GPL21145_pheno[, list(Sample_ID = title, Age = `age:ch1`)], by = "Sample_ID")
# subset(GPL21145_PhysAge, select = c("Sample_ID", "Age", setdiff(names(scores_all), "Sample_ID")))

### Output includes the following versions for the 8 surrogates and DNAmPhysAge:
# raw suffix indicates calculated scores of weighted CpGs
# _final suffix indicates the scores scaled in original values
# z prefix indicates standardized value (mean = 0, sd = 1)
# _years suffix indicates score is converted to age in years metric
