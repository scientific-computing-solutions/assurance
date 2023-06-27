context("negative binomial")

test_that("fit.nb matches MASS", {
  require(MASS)

  ## we'll use the quine dataset from MASS
  data <- with(
    quine,
    data.frame(
      days = Days,
      flag = as.numeric(Sex == "M"),
      duration = 1
    )
  )
  fit <- summary(glm.nb(days ~ offset(log(duration)) + flag,
    data = data, link = log
  ))
  my.fit <- with(
    quine,
    fit.nb(Days[Sex == "F"], Days[Sex == "M"], 1)
  )

  coefficients <- coef(fit)

  ## check that we match values of the regression model
  expect_equal(coefficients[, "Estimate"], my.fit$par[1:2],
    check.attributes = FALSE
  )

  ## and that we match the dispersion parameter
  expect_equal(fit$theta, 1 / my.fit$par[[3]], tol = 1e-6)

  ## and that we match standard errors
  expect_equal(coefficients[, "Std. Error"],
    sqrt(diag(my.fit$cov)),
    check.attributes = FALSE,
    tol = 1e-7
  )
})
