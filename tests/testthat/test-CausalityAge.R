test_that("calcCausalityAge functions work correctly", {
  skip_if_not(file.exists(test_path("fixtures", "CausalityAge", "GPL21145_matrix_CausalityAge.qs2")))
  GPL21145_matrix <- qs2::qs_read(
    test_path("fixtures", "CausalityAge", "GPL21145_matrix_CausalityAge.qs2")
  )
  GPL21145_CausAge_fixtures <- readRDS(
    test_path("fixtures", "CausalityAge", "GPL21145_CausAge.rds")
  )
  GPL21145_AdaptAge_fixtures <- readRDS(
    test_path("fixtures", "CausalityAge", "GPL21145_AdaptAge.rds")
  )
  GPL21145_DamAge_fixtures <- readRDS(
    test_path("fixtures", "CausalityAge", "GPL21145_DamAge.rds")
  )
  CausAge <- calcCausAge(GPL21145_matrix)
  AdaptAge <- calcAdaptAge(GPL21145_matrix)
  DamAge <- calcDamAge(GPL21145_matrix)
  expect_equal(CausAge$CausAge, GPL21145_CausAge_fixtures$CausAge)
  expect_equal(AdaptAge$AdaptAge, GPL21145_AdaptAge_fixtures$AdaptAge)
  expect_equal(DamAge$DamAge, GPL21145_DamAge_fixtures$DamAge)
})
