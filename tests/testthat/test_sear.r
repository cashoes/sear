library(sear)
context("Inputs")

test_that("Inputs are checked to contain meaningful gene names", {
  expect_warning(sear("ACTB"), "Submitted <10 genes. Enrichment may not be meaningful with so few genes as input.")
  expect_error(sear("actb"), "Make sure input are valid gene symbols.")
})
