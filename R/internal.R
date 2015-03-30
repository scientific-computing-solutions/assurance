# common utility functions used herein

#' @name internal
assert <- function(condition, msg="assertion failed") {
  # assert the condition to be true
  #
  # Args:
  #   condition : a boolean that must be TRUE
  #   msg       : a message to stop with
  if (!condition) {
    stop(msg)
  }
}

#' @name internal
is.wholenumber <- function(num) {
  # test to see if the argument is a whole number
  #
  # Args:
  #  num : the value(s) to test
  #
  # Returns:
  #   boolean
  is.integer(num) ||
  (is.numeric(num) && identical(num, round(num), num.eq=FALSE))
}

#' @name internal
is.positive <- function(num) {
  # test to see if the argument is positive
  #
  # Args:
  #   num : the value(s) to test
  #
  # Returns:
  #   single boolean
  is.numeric(num) && all(num >= 0)
}

#' @name internal
is.positive.scalar <- function(x) {
  is.numeric(x) && (length(x) == 1) && (x >= 0)
}

# test to see if the argument is a natural number (including zero)
# @param num the value(s) to test
# @return TRUE or FALSE (or NA)
#' @name internal
is.natural <- function(num) {
  is.wholenumber(num) && all(num >= 1)
}

#' @name internal
is.probability <- function(value) {
  is.numeric(value) && (min(value) >= 0) && (max(value) <= 1)
}


## Check the significance of this result
#' @name internal
check.significance <- function(alpha, expr, frame=parent.frame()) {
  if (is.na(alpha)) {
    TRUE
  } else {
    eval(expr, frame)
  }
}

check.hurdle <- function(values, possible.hurdle) {
  if (is.na(possible.hurdle)) {
    TRUE
  } else {
    values > possible.hurdle
  }
}

##' Fit a negative binomial regression model to two-arm count data.
##' @param comparatorEvents count of events for each patient on comparator
##' @param testEvents count of events for each patient on active
##' @param duration duration of trial
fit.nb <- function(comparatorEvents, testEvents, duration) {
  ## can get the treatment effect on rate exactly
  meanIntercept <- log(mean(comparatorEvents) / duration)
  meanTest      <- log(mean(testEvents) / duration) - meanIntercept

  ## wrap up everything for C++ code
  count <- c(comparatorEvents, testEvents)
  flag  <- c(rep(0, length(comparatorEvents)),
             rep(1, length(testEvents)))
  log.lambda <- meanIntercept + flag * meanTest

  ## get initial guess of dispersion parameter
  mc <- mean(count)
  eps.start <- (var(count) - mc) / mc^2
  coeffs <- c(max(0.01, eps.start))

  ## and now fit
  within(.Call(.call.fitNegBinData,
               coeffs, count, flag, log.lambda, duration),
         {
           par <- c(meanIntercept, meanTest, par)
         })

}
