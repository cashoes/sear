library(sear)
context("Inputs")

test_that("Inputs are checked to contain valid symbols", {
  expect_warning(sear("actb"))
})

test_that("Inputs are checked to contain sufficient number of valid symbols", {
  expect_warning(sear("ACTB"))
  expect_error(sear(c("ACTB", "foo", "bar", "baz", "bat")))
})
