context("utility functions")


test_that("assert works", {
  expect_error(assert(FALSE, "test"))
  expect_true(is.null(assert(TRUE, "test")))
})


test_that("is.wholenumber works", {
  expect_true(is.wholenumber(1))
  expect_false(is.wholenumber(1.1))
  expect_false(is.wholenumber("lemon"))
  expect_true(is.wholenumber(1:10))
})


test_that("is.positive works", {
  expect_true(is.positive(0.0))
  expect_false(is.positive(-1.0))
  expect_false(is.positive("custard"))
  expect_true(is.positive(10.0))
})


test_that("is.natural works", {
  expect_true(is.natural(1))
  expect_false(is.natural(0))
  expect_false(is.natural(1.1))
  expect_false(is.natural(-1))
})


test_that("study.size works", {
  expect_error(study.size(total.size = -1))
  expect_error(study.size(total.size = 0))
  expect_error(study.size(total.size = 1))
  expect_error(study.size())
  expect_error(study.size(grp2.size = 2))
  expect_error(study.size(grp1.size = -1, grp2.size = -2))
  fred <- study.size(grp1.size = 10, grp2.size = 20)
  expect_equal(fred@grp1, 10)
  expect_equal(fred@grp2, 20)
  fred <- study.size(total.size = 20)
  expect_equal(fred@grp1, 10)
  expect_equal(fred@grp2, 10)
})


test_that("pts.rec works in simple cases", {
  expect_equal(
    pts.rec(0, 2, c(0.1, 0.2)),
    0.9 * 0.8
  )
  expect_equal(
    pts.rec(1, 2, c(0.1, 0.2)),
    0.1 * 0.8 + 0.9 * 0.2
  )
  expect_equal(pts.rec(2, 2, c(0.1, 0.2)), 0.1 * 0.2)
})

test_that("pts.rec recurrence works", {
  set.seed(20150323)
  probs <- runif(5)

  three <- pts.rec(3, 5, probs)

  for (ind in seq_along(probs)) {
    dropped <- probs[-ind]
    rec <- (1 - probs[[ind]]) * pts.rec(3, 4, dropped) +
      probs[[ind]] * pts.rec(2, 4, dropped)
    expect_equal(three, rec)
  }
})
