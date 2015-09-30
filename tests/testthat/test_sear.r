library(sear)
context("Inputs")

test_that("Inputs are checked to contain meaningful gene names", {
  expect_warning(sear("ACTB", type = "mrna"), "You submitted <10 valid symbols. Results may not be meaningful with so few genes.")
  expect_error(sear("actb", type = "mrna"), "You selected type = %s, but many of your features are not recognized.")
})
