context("single arm survival")

test_that("known value works", {
  ## this is a case that Andy Stone did in a spreadsheet
  ## so it's an independent implementation
  initial <- new.singleArm(18, 48, 1, 1)
  rate <- new.oneArmResponseRate(77, 0.29)
  result <- assurance(initial, rate)
  ## and here's the result.
  expect_equal(result, 0.84057583)
})
