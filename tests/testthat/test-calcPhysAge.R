test_that("calcPhysAge", {
  skip_if_not(file.exists(test_path("fixtures", "PhysAge", "GPL21145_PhysAge.rds")))
  GPL21145_fixtures <- readRDS(test_path("fixtures", "PhysAge", "GPL21145_PhysAge.rds"))

  GPL21145_matrix <- qs2::qs_read(test_path("fixtures", "PhysAge", "GPL21145_matrix_PhysAge.qs2"))
  PhysAge <- calcPhysAge(GPL21145_matrix)
  PhysAge <- PhysAge[match(GPL21145_fixtures$Sample_ID, PhysAge$Sample_ID), names(GPL21145_fixtures)]

  row.names(GPL21145_fixtures) <- NULL
  row.names(PhysAge) <- NULL

  expect_equal(PhysAge, GPL21145_fixtures)
})

