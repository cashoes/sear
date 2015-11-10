library(sear)
context("Inputs")

test_that("Inputs are checked to contain meaningful gene names", {
  expect_warning(sear("ACTB"), "You submitted <10 valid symbols. Results may not be meaningful with so few inputs.")
  expect_error(sear("actb"))
  expect_error(sear(c("ACTB", "foo", "bar", "baz", "bat")))
})
