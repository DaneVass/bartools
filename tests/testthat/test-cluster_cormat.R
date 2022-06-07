test_that("cluster cormat works", {
  test.mat <- readRDS("../test.mat.rds")
  test.mat.clust <- readRDS("../test.mat.clust.rds")
  expect_equal(cluster_cormat(test.mat), test.mat.clust)
})
