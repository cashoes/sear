library(sear)
context("Inputs")

test_that("Inputs are checked to contain valid symbols", {
  expect_warning(sear("ACTB"))
  expect_error(sear("actb"))
  expect_error(sear(c("ACTB", "foo", "bar", "baz", "bat")))
})

test_that("Inputs are checked to contain the correct type of input", {
  expect_warning(sear("ACTB"))
  expect_error(sear("actb"))
  expect_error(sear(c("ACTB", "foo", "bar", "baz", "bat")))
})

test_that("Inputs are checked to contain sufficient number of valid symbols", {
  expect_warning(sear("ACTB"))
  expect_error(sear("actb"))
  expect_error(sear(c("ACTB", "foo", "bar", "baz", "bat")))
})
