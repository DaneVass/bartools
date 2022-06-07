test_that("get_upper_tri works", {
  test.mat <- readRDS("../test.mat.rds")
  test.mat.upper <- readRDS("../test.mat.upper.rds")
  expect_equal(get_upper_tri(test.mat), test.mat.upper)
})
