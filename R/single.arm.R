#' @include assurance.R
NULL

##' Assurance for single arm response rate
##'
##' Assurance can be calculated for response rate
##'
##' Arbitrary beta distributions can be used as prior.
##'
##' @aliases length,singleArmProbabilities-method
##' samplePosterior,singleArmInitialData,numeric-method
##' treatmentEffect,singleArmProbabilities-method
##' sampleLater,singleArmProbabilities,oneArmResponseRate-method
##' testLater,singleArmLater-method
##' sampleAndTestLater,singleArmProbabilities,oneArmResponseRate-method
##' assurance,singleArmInitialData,oneArmResponseRate,missing-method
##' @name single-arm
##' @examples
##' initial.data <- new.singleArm(15, 30)
##' later <- new.oneArmResponseRate(100, 0.43)
##' # do it numerically
##' assurance(initial.data, later, 1000000)
##' # do it analytically
##' assurance(initial.data, later)
NULL

setClass("singleArmInitialData",
  representation(
    num.success = "numeric",
    total.num = "numeric",
    shape1 = "numeric",
    shape2 = "numeric"
  ),
  validity = function(object) {
    if (!is.positive.scalar(object@num.success)) {
      "num.success must be a positive scalar"
    } else if (!is.positive.scalar(object@total.num)) {
      "total.num must be a positive scalar"
    } else if (object@num.success > object@total.num) {
      "num.success must be less than or equal to total.num"
    } else {
      TRUE
    }
  }
)

##' @export new.singleArm
new.singleArm <- function(num.success, total.num, shape1 = 0.5, shape2 = 0.5) {
  shape1 <- num.success + shape1
  shape2 <- total.num - num.success + shape2
  new("singleArmInitialData",
    num.success = num.success,
    total.num = total.num,
    shape1 = shape1,
    shape2 = shape2
  )
}


setClass("singleArmProbabilities",
  representation(cure.probability = "numeric"),
  validity = function(object) {
    if (!is.probability(object@cure.probability)) {
      "cure.probability must be in [0, 1]"
    } else {
      TRUE
    }
  }
)


setMethod(
  "length",
  signature(x = "singleArmProbabilities"),
  function(x) {
    length(x@cure.probability)
  }
)


setMethod(
  "samplePosterior",
  signature(
    earlyStudy = "singleArmInitialData",
    nsim = "numeric"
  ),
  function(earlyStudy, nsim) {
    prob <- rbeta(nsim, earlyStudy@shape1, earlyStudy@shape2)
    new("singleArmProbabilities", cure.probability = prob)
  }
)

setMethod(
  "treatmentEffect",
  signature(posteriorSample = "singleArmProbabilities"),
  function(posteriorSample) {
    posteriorSample@cure.probability
  }
)

setClass(
  "singleArmLater",
  representation(
    success.count = "numeric",
    study.defn = "oneArmResponseRate"
  )
)


setMethod(
  "sampleLater",
  signature(
    posteriorSample = "singleArmProbabilities",
    laterStudy = "oneArmResponseRate"
  ),
  function(posteriorSample, laterStudy) {
    nsim <- length(posteriorSample)
    success.count <- rbinom(
      nsim,
      laterStudy@size@size,
      posteriorSample@cure.probability
    )
    new("singleArmLater",
      success.count = success.count,
      study.defn = laterStudy
    )
  }
)


setMethod(
  "testLater",
  signature(laterSample = "singleArmLater"),
  function(laterSample) {
    rate <- laterSample@study.defn@rate
    size <- laterSample@study.defn@size@size
    wanted <- ceiling(rate * size)
    laterSample@success.count >= wanted
  }
)


setMethod(
  "sampleAndTestLater",
  signature(
    posteriorSample = "singleArmProbabilities",
    laterStudy = "oneArmResponseRate"
  ),
  function(posteriorSample, laterStudy) {
    rate <- laterStudy@rate
    size <- laterStudy@size@size
    wanted <- ceiling(rate * size)
    pbinom(size - wanted, size,
      1 - posteriorSample@cure.probability,
      lower.tail = TRUE
    )
  }
)

setMethod(
  "assurance",
  signature(
    earlyStudy = "singleArmInitialData",
    laterStudy = "oneArmResponseRate",
    nsim = "missing"
  ),
  function(earlyStudy, laterStudy, nsim) {
    rate <- laterStudy@rate
    size <- laterStudy@size@size
    wanted <- ceiling(rate * size)
    k <- seq(wanted, size)
    bin.coeffs <- lchoose(size, k)
    top <- lbeta(
      earlyStudy@shape1 + k,
      earlyStudy@shape2 + size - k
    )
    bottom <- lbeta(earlyStudy@shape1, earlyStudy@shape2)
    sum(exp(bin.coeffs + top - bottom))
  }
)
