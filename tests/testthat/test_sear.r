library(sear)
context("Inputs")

test_that("Inputs are checked to contain valid symbols", {
  expect_error(sear("actb"))
})

test_that("Inputs are checked to contain sufficient number of valid symbols", {
  expect_warning(sear(c("ACTB", "foo", "bar", "baz", "bat")))
})
