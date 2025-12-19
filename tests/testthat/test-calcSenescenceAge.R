test_that("calcSenescenceAge functions work correctly", {
  skip_if_not(file.exists(test_path("fixtures", "SenescenceAge", "GPL21145_matrix_SenescenceAge.qs2")))
  GPL21145_matrix <- qs2::qs_read(
    test_path("fixtures", "SenescenceAge", "GPL21145_matrix_SenescenceAge.qs2")
  )

  GPL21145_SenChronoAge_fixtures <- readRDS(
    test_path("fixtures", "SenescenceAge", "GPL21145_SenChronoAge.rds")
  )
  GPL21145_SenCultureAge_fixtures <- readRDS(
    test_path("fixtures", "SenescenceAge", "GPL21145_SenCultureAge.rds")
  )
  GPL21145_SenMortalityAge_fixtures <- readRDS(
    test_path("fixtures", "SenescenceAge", "GPL21145_SenMortalityAge.rds")
  )

  SenChronoAge <- calcSenChronoAge(GPL21145_matrix)
  SenCultureAge <- calcSenCultureAge(GPL21145_matrix)
  SenMortalityAge <- calcSenMortalityAge(GPL21145_matrix)

  expect_equal(SenChronoAge, GPL21145_SenChronoAge_fixtures)
  expect_equal(SenCultureAge, GPL21145_SenCultureAge_fixtures)
  expect_equal(SenMortalityAge, GPL21145_SenMortalityAge_fixtures)
})
