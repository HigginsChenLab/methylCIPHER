library(data.table)

PhysAge_CpGs <- fread("data-raw/PhysAge/CpG lists and weights/HRS_cpg_training_means.txt", header = TRUE)
use_data(PhysAge_CpGs, overwrite = TRUE)

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

PhysAge_data <- lapply(weight_files, \(x) fread(x, header = FALSE, col.names = c("CpG", "Beta")))
use_data(PhysAge_data, overwrite = TRUE)
