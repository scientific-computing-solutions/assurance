context("survival")

test_that("object creation works", {
  expect_error(new.survival(0.3, -1, 10))
  fred <- new.survival(0.3, 100, 5000)
  expect_equal(fred@size@grp1, 100)
  expect_equal(fred@size@grp2, 5000)
  expect_equal(fred@log.hazard.ratio, log(0.3))
})


test_that("sampling works", {
  fred <- new.survival(0.3, 100, 5000)
  nsim <- 100000
  sample <- samplePosterior(fred, nsim)
  expect_equal(length(sample), nsim)
  expect_equal(mean(sample@log.hazard.ratio), log(0.3), tolerance = 1e-2)
})
