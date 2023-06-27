#' @include assurance.R
NULL

## assurance for (overdispersed) count data

setClass("poissonOverDispersedArm",
  representation(
    count = "numeric",
    total.exposure = "numeric",
    overdispersion = "numeric"
  ),
  validity = function(object) {
    if (!is.positive(object@total.exposure)) {
      "total.exposure must be positive"
    } else if (!is.positive(object@count)) {
      "count must be positive"
    } else if (!is.positive(object@overdispersion)) {
      "overdispersion parameter must be positive"
    } else {
      TRUE
    }
  }
)

setMethod(
  "show",
  signature(object = "poissonOverDispersedArm"),
  function(object) {
    cat("Event count   : ", object@count, "\n", sep = "")
    cat("Time at risk  : ", object@total.exposure, "\n", sep = "")
    cat("Overdispersion: ", object@overdispersion, "\n", sep = "")
    invisible(NULL)
  }
)

setClass(
  "poissonInitialData",
  representation(
    grp1.data = "poissonOverDispersedArm",
    grp2.data = "poissonOverDispersedArm"
  )
)

setMethod(
  "show",
  signature(object = "poissonInitialData"),
  function(object) {
    cat("Overdispersed Poisson\n\n")
    cat("Group 1:\n")
    show(object@grp1.data)
    cat("\n")
    cat("Group 2:\n")
    show(object@grp2.data)
    invisible(NULL)
  }
)

## here's the factory function
##' @export new.overdispersedPoisson
new.overdispersedPoisson <-
  function(grp1.overdispersion, grp1.count, grp1.time.at.risk,
           grp2.overdispersion, grp2.count, grp2.time.at.risk) {
    new("poissonInitialData",
      grp1.data = new("poissonOverDispersedArm",
        count = grp1.count,
        total.exposure = grp1.time.at.risk,
        overdispersion = grp1.overdispersion
      ),
      grp2.data = new("poissonOverDispersedArm",
        count = grp2.count,
        total.exposure = grp2.time.at.risk,
        overdispersion = grp2.overdispersion
      )
    )
  }

samplePoissonRate <- function(obj, nsim) {
  ## which had better be a poissonOverDispersedArm
  num.events <- obj@count
  total.time.at.risk <- obj@total.exposure
  overdispersion <- obj@overdispersion
  alpha <- num.events / overdispersion
  beta <- total.time.at.risk / overdispersion
  rgamma(nsim, shape = alpha, rate = beta)
}

setClass("poissonRates",
  representation(
    grp1.rate = "numeric",
    grp2.rate = "numeric"
  ),
  validity = function(object) {
    r1 <- object@grp1.rate
    r2 <- object@grp2.rate
    if (length(r1) != length(r2)) {
      "rates must have same length"
    } else if (any(r1 <= 0) || any(r2 <= 0)) {
      "rates must be strictly positive"
    } else {
      TRUE
    }
  }
)

setMethod(
  "length",
  signature(x = "poissonRates"),
  function(x) {
    length(x@grp1.rate)
  }
)

setMethod(
  "treatmentEffect",
  signature(posteriorSample = "poissonRates"),
  function(posteriorSample) {
    ## the treatment effect is the relative rate reduction
    r1 <- posteriorSample@grp1.rate
    r2 <- posteriorSample@grp2.rate
    1 - (r1 / r2)
  }
)

setMethod(
  "samplePosterior",
  signature(earlyStudy = "poissonInitialData", nsim = "numeric"),
  function(earlyStudy, nsim) {
    grp1.rate <- samplePoissonRate(earlyStudy@grp1.data, nsim)
    grp2.rate <- samplePoissonRate(earlyStudy@grp2.data, nsim)
    new("poissonRates",
      grp1.rate = grp1.rate,
      grp2.rate = grp2.rate
    )
  }
)


setClass("overdispersedPoissonStudy",
  contains = "twoArm",
  representation(
    grp1.overdispersion = "numeric",
    grp2.overdispersion = "numeric",
    duration = "numeric"
  ),
  validity = function(object) {
    if (!is.positive(object@grp1.overdispersion)) {
      "grp1.overdispersion must be positive"
    } else if (!is.positive(object@grp2.overdispersion)) {
      "grp2.overdispersion must be positive"
    } else if (!is.positive(object@duration)) {
      "duration must be positive"
    } else {
      TRUE
    }
  }
)


## here's the factory function
##' @export new.overdispersedPoissonStudy
new.overdispersedPoissonStudy <-
  function(duration, grp1.size, grp2.size,
           grp1.overdispersion, grp2.overdispersion,
           significance = NA,
           hurdle = NA) {
    new("overdispersedPoissonStudy",
      duration = duration,
      grp1.overdispersion = grp1.overdispersion,
      grp2.overdispersion = grp2.overdispersion,
      significance = as.numeric(significance),
      hurdle = as.numeric(hurdle),
      size = study.size(grp1.size = grp1.size, grp2.size = grp2.size)
    )
  }

## get total time-at-risk in both groups
durations <- function(study) {
  size <- study@size
  duration <- study@duration
  list(
    grp1 = size@grp1 * duration,
    grp2 = size@grp2 * duration
  )
}

setClass("overdispersedPoissonLater",
  representation(
    grp1.count = "numeric",
    grp2.count = "numeric",
    study.defn = "overdispersedPoissonStudy"
  ),
  validity = function(object) {
    if (length(object@grp1.count) != length(object@grp2.count)) {
      "counts must be same length"
    } else if (any(object@grp1.count < 0) ||
      any(object@grp2.count < 0)) {
      "counts must be nonnegative"
    } else {
      TRUE
    }
  }
)

setMethod(
  "sampleLater",
  signature(
    posteriorSample = "poissonRates",
    laterStudy = "overdispersedPoissonStudy"
  ),
  function(posteriorSample, laterStudy) {
    duration <- durations(laterStudy)
    nsim <- length(posteriorSample)
    grp1.count <- rpois(
      nsim,
      posteriorSample@grp1.rate * duration$grp1
    )
    grp2.count <- rpois(
      nsim,
      posteriorSample@grp2.rate * duration$grp2
    )
    new("overdispersedPoissonLater",
      grp1.count = grp1.count,
      grp2.count = grp2.count,
      study.defn = laterStudy
    )
  }
)


setMethod(
  "testLater",
  signature(laterSample = "overdispersedPoissonLater"),
  function(laterSample) {
    study.def <- laterSample@study.defn
    alpha <- study.def@significance
    duration <- durations(study.def)
    hurdle <- study.def@hurdle
    if (any(laterSample@grp1.count == 0) ||
      any(laterSample@grp2.count == 0)) {
      stop("Zero counts: inappropriate test")
    }
    log.rate.1 <- log(laterSample@grp1.count) - log(duration$grp1)
    log.rate.2 <- log(laterSample@grp2.count) - log(duration$grp2)
    effect <- log.rate.1 - log.rate.2
    check.significance(
      alpha,
      {
        z <- qnorm(0.5 * alpha, lower.tail = FALSE)
        se <-
          sqrt(study.def@grp1.overdispersion / laterSample@grp1.count +
            study.def@grp2.overdispersion / laterSample@grp2.count)
        abs(effect) > se * z
      }
    ) &
      check.hurdle(-expm1(effect), hurdle)
  }
)


setClass("negativeBinomialStudy",
  contains = "twoArm",
  representation(
    overdispersion = "numeric",
    duration = "numeric"
  ),
  validity = function(object) {
    if (!is.positive(object@overdispersion)) {
      "overdispersion must be positive"
    } else if (!is.positive(object@duration)) {
      "duration must be positive"
    } else {
      TRUE
    }
  }
)

##' @export new.negativeBinomialStudy
new.negativeBinomialStudy <-
  function(duration, grp1.size, grp2.size, overdispersion,
           significance = NA, hurdle = NA) {
    new("negativeBinomialStudy",
      duration = duration,
      overdispersion = overdispersion,
      significance = as.numeric(significance),
      hurdle = as.numeric(hurdle),
      size = study.size(grp1.size = grp1.size, grp2.size = grp2.size)
    )
  }

setClass(
  "negativeBinomialLater",
  representation(
    rates = "poissonRates",
    study.defn = "negativeBinomialStudy"
  )
)

## we don't have summary stats here so we're going to have to simulate
## the actual trial
setMethod(
  "sampleLater",
  signature(
    posteriorSample = "poissonRates",
    laterStudy = "negativeBinomialStudy"
  ),
  function(posteriorSample, laterStudy) {
    new("negativeBinomialLater",
      rates = posteriorSample,
      study.defn = laterStudy
    )
  }
)

setMethod(
  "testLater",
  signature(laterSample = "negativeBinomialLater"),
  function(laterSample) {
    study.def <- laterSample@study.defn
    overdispersion <- study.def@overdispersion
    rates <- rbind(
      laterSample@rates@grp1.rate,
      laterSample@rates@grp2.rate
    )
    flag <- c(
      rep(0, study.def@size@grp1),
      rep(1, study.def@size@grp2)
    )
    apply(
      rates, 2,
      function(column, dispersion, flag,
               grp1.size, grp2.size, duration,
               alpha, hurdle) {
        fuzz.1 <- rgamma(grp1.size,
          shape = dispersion,
          rate = dispersion
        )
        fuzz.2 <- rgamma(grp2.size,
          shape = dispersion,
          rate = dispersion
        )
        count.1 <- rpois(
          grp1.size,
          fuzz.1 * column[[1]] * duration
        )
        count.2 <- rpois(
          grp2.size,
          fuzz.2 * column[[2]] * duration
        )
        fit <- fit.nb(count.2, count.1, duration)
        effect <- fit$par[[2]]
        check.significance(
          alpha,
          {
            z <- qnorm(alpha, lower.tail = FALSE)
            se <- sqrt(fit$cov[2, 2])
            abs(effect) > z * se
          }
        ) &
          check.hurdle(-expm1(effect), hurdle)
      },
      1 / overdispersion,
      flag,
      study.def@size@grp1,
      study.def@size@grp2,
      study.def@duration,
      0.5 * study.def@significance,
      study.def@hurdle
    )
  }
)
