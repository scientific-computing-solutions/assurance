context("gaussian endpoints")

test_that("rsichisq works", {
  ## test scaled inverse chi^2
  ## according to Gelman p575
  ## we check the expectation matches the Gelman expectation
  df <- 10
  scale <- 5
  nsim <- 10000000
  samples <- rsichisq(nsim, df, scale)
  expect_equal(length(samples), nsim)
  expectation <- df / (df - 2) * scale
  expect_equal(mean(samples), expectation, tolerance = 1e-2)
})


test_that("new.gaussian works", {
  expect_error(new.gaussian(10, -1, -1, sd = 10))
  expect_error(new.gaussian(10, 1.1, 1.1, sd = 10))
  expect_error(new.gaussian(10, "custard", 1.1))
  got <- new.gaussian(10, 5, 10, sd = 15)
  expect_equal(got@size@grp1, 5)
  expect_equal(got@size@grp2, 10)
  expect_equal(got@delta.mu, 10)
  expect_equal(got@sigma2, 225)
})


test_that("common sd sampling works", {
  got <- new.gaussian(3, 10, 10, sd = 8)
  nsim <- 1000000
  sample <- samplePosterior(got, nsim)
  delta.mu <- sample@delta.mu
  expect_equal(got@delta.mu, 3)
  expect_equal(got@sigma2, 64)
  expect_equal(length(delta.mu), nsim)
  expect_equal(length(sample@sigma2), nsim)
  expect_equal(mean(delta.mu), 3, tolerance = 1e-2)
})


test_that("different sd sampling works", {
  got <- new.gaussian(3, 20, 50, sd1 = 4, sd2 = 10)
  expect_equal(got@delta.mu, 3)
  expect_equal(got@grp1.sigma2, 16)
  expect_equal(got@grp2.sigma2, 100)
  expect_equal(got@size@grp1, 20)
  expect_equal(got@size@grp2, 50)
  nsim <- 100000
  sample <- samplePosterior(got, nsim)
  expect_equal(length(sample@delta.mu), nsim)
  expect_equal(length(sample@grp1.sigma2), nsim)
  expect_equal(length(sample@grp2.sigma2), nsim)
})
