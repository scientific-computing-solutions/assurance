context("binomial endpoints")

test_that("new.binomial throws errors", {
  expect_error(new.binomial(x1=-1, m1=10, x2=0, m2=10))
  expect_error(new.binomial(x1=11, m1=10, x2=0, m2=10))
  expect_error(new.binomial(x2=-1, m2=10, x1=0, m1=10))
  expect_error(new.binomial(x2=11, m2=10, x1=0, m1=10))
})

test_that("sample size is right", {
  fred <- new.binomial(x1=5, m1=100, x2=10, m2=100)
  nsim <- 666
  sample <- samplePosterior(fred, nsim)
  expect_equal(length(sample), nsim)
})
