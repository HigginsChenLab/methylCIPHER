# SystemsAge_data ----
## This file convert the `SystemsAge_data.RData` into qs2 trimmed down to save space.
## There's actually not a lot to trim because most of the mass is in the rotation
## object which is a dense matrix. But saving to qs2 saved loading time by 5-6 times
## on a multithreaded system.

load_all()

system.time(
  load(
    paste0(
      get_methylCIPHER_default_path(),
      "/",
      "SystemsAge_data.RData"
    )
  )
)

#  user  system elapsed
# 11.37    1.06   18.08

# names(DNAmPCA)
# DNAmPCA |> lapply(lobstr::obj_size)
all.equal(CpGs, names(imputeMissingCpGs))

SystemsAge_data <- list(
  CpGs = CpGs,
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

lapply(SystemsAge_data, lobstr::obj_size)

# generate the intended hash
rlang::hash(SystemsAge_data)

qs2::qs_save(
  SystemsAge_data,
  file.path(get_methylCIPHER_path(), "SystemsAge_data.qs2")
)

system.time(
  SystemsAge_data <- qs2::qs_read(
    file.path(get_methylCIPHER_path(), "SystemsAge_data.qs2")
  )
)

# user  system elapsed
# 2.06    0.67    2.80
