test_that("calcSystemsAge works", {
  skip("Manual Test Only")
  expect_no_error(
    suppressWarnings(calcSystemsAge(DNAm = exampleBetas, pheno = examplePheno))
  )
})
