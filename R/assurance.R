#' @include internal.R studies.R
NULL

setGeneric("samplePosterior",
  def = function(earlyStudy, nsim) {
    standardGeneric("samplePosterior")
  }
)

setGeneric("sampleInitialTrial",
  def = function(earlyStudy, nsim) {
    standardGeneric("sampleInitialTrial")
  }
)


setGeneric("sampleLater",
  def = function(posteriorSample, laterStudy) {
    standardGeneric("sampleLater")
  }
)


setGeneric("testLater",
  def = function(laterSample) standardGeneric("testLater")
)

setGeneric("sampleAndTestLater",
  def = function(posteriorSample, laterStudy) {
    testLater(sampleLater(posteriorSample, laterStudy))
  }
)


setGeneric("treatmentEffect",
  def = function(posteriorSample) {
    standardGeneric("treatmentEffect")
  }
)

setGeneric("assurance",
  useAsDefault = function(earlyStudy, laterStudy, nsim) {
    posterior <- samplePosterior(earlyStudy, nsim)
    assert(length(posterior) == nsim, "not enough samples")
    mean(sampleAndTestLater(posterior, laterStudy))
  }
)



setGeneric("assurance.multiple.endpoints",
  useAsDefault = function(earlyStudy, laterStudy, nsim) {
    posterior <- samplePosterior(earlyStudy, nsim)
    # assert( length( posterior ) == nsim, "not enough samples" )
    sampleAndTestLater(posterior, laterStudy)
  }
)


# here are methods that handle lists of later studies
# this is cunning code
# so be warned


setClass(
  "generalLaterStudy",
  representation(
    required.success = "numeric",
    studies = "list",
    at.least = "logical"
  )
)


##' @export new.generalLaterStudy
new.generalLaterStudy <- function(required.success, studies, at.least = F) {
  stopifnot(required.success >= 0)
  stopifnot(length(studies) > 0)
  new("generalLaterStudy", required.success = required.success, studies = studies, at.least = at.least)
}

setClass(
  "generalLaterSample",
  representation(
    required.success = "numeric",
    results = "list",
    at.least = "logical"
  )
)

setMethod(
  "sampleLater",
  signature(
    posteriorSample = "ANY",
    laterStudy = "generalLaterStudy"
  ),
  function(posteriorSample, laterStudy) {
    new("generalLaterSample",
      required.success = laterStudy@required.success,
      results = lapply(
        laterStudy@studies,
        function(x) {
          sampleLater(posteriorSample, x)
        }
      ),
      at.least = laterStudy@at.least
    )
  }
)



setMethod(
  "sampleLater",
  signature(
    posteriorSample = "ANY",
    laterStudy = "list"
  ),
  function(posteriorSample, laterStudy) {
    lapply(
      laterStudy,
      function(x) {
        sampleLater(posteriorSample, x)
      }
    )
  }
)

## ! Function for calculating probability of having k out of n technical successes
## ! @param k integer
## ! @param n integer
## ! @return probability of having k out of n technical successes
pts.rec <- function(k, n, p) {
  stopifnot(k <= n)
  stopifnot(n > 0)
  stopifnot(length(p) >= n)

  q <- 1 - p

  # Special case
  if (k == 0) {
    return(prod(q[1:n]))
  }

  # Base case
  if (k == 1 & n == 1) {
    return(p[1])
  }

  # Case 1: k = n
  if (k == n) {
    return(prod(p[1:n]))
  }

  # Case 2: k=1, n>1
  if (k == 1) {
    return((1 - p[n]) * pts.rec(1, n - 1, p) + p[n] * prod(q[1:n - 1]))
  }

  # Case 3: k<n and k>1
  if (k < n) {
    return((1 - p[n]) * pts.rec(k, n - 1, p) + p[n] * pts.rec(k - 1, n - 1, p))
  }

  warning("Something is not right")
  return(NA)
}


setMethod(
  "testLater",
  signature(laterSample = "generalLaterSample"),
  function(laterSample) {
    results <- sapply(laterSample@results, testLater)
    required.success <- laterSample@required.success

    res <- 0
    if (laterSample@at.least == F) {
      # Calculate k out of n
      res <- apply(results, 1, pts.rec, k = required.success, n = ncol(results))
    } else {
      # Calculate at least k out of n
      for (i in required.success:ncol(results)) {
        res <- res + apply(results, 1, pts.rec, k = i, n = ncol(results))
      }
    }
    return(res)
  }
)



setMethod(
  "testLater",
  signature(laterSample = "list"),
  function(laterSample) {
    # so, get all the test results
    test.results <- sapply(laterSample, testLater)
    # and for a given posterior, we want all to succeed
    apply(test.results, 1, prod)
  }
)


setClass("twoArm",
  representation(
    size = "twoWayStudySize",
    significance = "numeric",
    hurdle = "numeric",
    test.me = "character",
    m1 = "numeric",
    m2 = "numeric"
  ),
  validity = function(object) {
    signif <- object@significance
    hurdle <- object@hurdle
    if (length(signif) != 1) {
      "significance must be a scalar"
    } else if (!(is.na(signif) || is.probability(signif))) {
      "significance must be a probability"
    } else if (is.na(signif) && is.na(hurdle)) {
      "must give significance or hurdle or both"
    } else {
      TRUE
    }
  }
)


##' @export new.twoArm
new.twoArm <- function(size, significance = NULL, hurdle = NULL,
                       test.me = "Uncorrected.Significant", m1 = NULL, m2 = NULL) {
  if (is.null(significance)) {
    significance <- as.numeric(NA)
  }
  if (is.null(hurdle)) {
    hurdle <- as.numeric(NA)
  }
  if (is.null(test.me)) {
    test.me <- as.numeric(NA)
  }
  if (is.null(m1)) {
    m1 <- as.numeric(NA)
  }
  if (is.null(m2)) {
    m2 <- as.numeric(NA)
  }
  new("twoArm", size = size, significance = significance, hurdle = hurdle, test.me = test.me, m1 = m1, m2 = m2)
}



setClass("oneArmResponseRate",
  representation(
    size = "oneWayStudySize",
    rate = "numeric"
  ),
  validity = function(object) {
    rate <- object@rate
    if (length(rate) != 1) {
      "significance must be a scalar"
    } else if (!is.probability(rate)) {
      "significance must be a probability"
    } else {
      TRUE
    }
  }
)


##' @export new.oneArmResponseRate
new.oneArmResponseRate <- function(size, rate) {
  new("oneArmResponseRate",
    size = new("oneWayStudySize", size = size),
    rate = rate
  )
}
