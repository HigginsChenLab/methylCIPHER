test_that("calcSystemsAge works", {
  expect_no_error(
    suppressWarnings(calcSystemsAge(DNAm = exampleBetas, pheno = examplePheno))
  )
})
