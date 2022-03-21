test_that("sort", {

  y <- runif(30)
  x <- matrix(runif(30 * 30), 30, 30)

  expect_equal(nrow(p.sort(x, y)), 30)
  expect_equal(ncol(p.sort(x, y)), 30)
})
