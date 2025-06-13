test_that("calcPCClocks works", {
  expect_no_error(
    suppressWarnings(calcPCClocks(DNAm = exampleBetas, pheno = examplePheno))
  )
})
