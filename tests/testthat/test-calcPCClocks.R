test_that("calcPCClocks works", {
  skip("Manual Test Only")
  expect_no_error(
    suppressWarnings(calcPCClocks(DNAm = exampleBetas, pheno = examplePheno))
  )
})
