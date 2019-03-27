
context("model_pod output")
library(testthat)
library(modelpod)


load("../testdata/testdata.RData")
result <- model_pod(xcol = "dose_uM", ycol = "GFP_int", dataset = test_data, respLev = 0.1, groups = "treatment")


test_that("model_pod returns dataframe", {
  expect_equal(class(result), "data.frame")
  expect_equal(nrow(result), 30)

})


test_that("Estimates are correct", {
  expect_equal((result[3,2] - 185.7214) < 0.001, TRUE )
  expect_equal(is.na(result[19,2]),  TRUE)

})


