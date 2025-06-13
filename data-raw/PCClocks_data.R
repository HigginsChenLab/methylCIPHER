# PCClocks_data ----
## This file convert the `CalcAllPCClocks.RData` into qs2 trimmed down to save space.
## There's actually not a lot to trim because most of the mass is in the rotation
## object which is a dense matrix. But saving to qs2 saved loading time by 5-6 times
## on a multithreaded system.

load_all()

system.time(
  load(
    paste0(
      get_methylCIPHER_default_path(),
      "/",
      "CalcAllPCClocks.RData"
    )
  )
)

#   user  system elapsed
#   9.88    0.44   11.48

all.equal(CpGs, names(imputeMissingCpGs))

PCClocks_data <- list(
  CalcPCDNAmTL = CalcPCDNAmTL,
  CalcPCPhenoAge = CalcPCPhenoAge,
  CalcPCGrimAge = CalcPCGrimAge,
  CalcPCHannum = CalcPCHannum,
  imputeMissingCpGs = imputeMissingCpGs,
  CalcPCHorvath1 = CalcPCHorvath1,
  CalcPCHorvath2 = CalcPCHorvath2
)

lapply(PCClocks_data, lobstr::obj_size)

# generate the intended hash
rlang::hash(PCClocks_data)

qs2::qs_save(
  PCClocks_data,
  file.path(get_methylCIPHER_path(), "PCClocks_data.qs2")
)

system.time(
  PCClocks_data <- qs2::qs_read(
    file.path(get_methylCIPHER_path(), "PCClocks_data.qs2")
  )
)

# user  system elapsed
# 1.28    0.56    2.28
