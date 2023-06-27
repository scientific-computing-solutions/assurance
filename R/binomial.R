#' @include assurance.R
NULL

setClass("binomialInitialData",
  representation(
    grp1.count = "numeric",
    grp2.count = "numeric",
    size = "twoWayStudySize"
  ),
  validity = function(object) {
    grp1.count <- object@grp1.count
    grp2.count <- object@grp2.count
    size <- object@size
    if (!(is.wholenumber(grp1.count) &&
      is.positive(grp1.count))) {
      "grp1.count must be a positive whole number"
    } else if (!(is.wholenumber(grp2.count) &&
      is.positive(grp2.count))) {
      "grp2.count must be a positive whole number"
    } else if (grp1.count > size@grp1) {
      "grp1.count must be <= grp1.size"
    } else if (grp2.count > size@grp2) {
      "grp2.count must be <= grp2.size"
    } else {
      TRUE
    }
  }
)


setClass("binomialProbabilities",
  representation(
    grp1.prob = "numeric",
    grp2.prob = "numeric"
  ),
  validity = function(object) {
    p1 <- object@grp1.prob
    p2 <- object@grp2.prob
    if (length(p1) != length(p2)) {
      "grp1.prob and grp2.prob must have same lengths"
    } else if (!is.probability(p1)) {
      "grp1.prob must be in [0, 1]"
    } else if (!is.probability(p2)) {
      "grp2.prob must be in [0, 1]"
    } else {
      TRUE
    }
  }
)


setMethod(
  "length",
  signature(x = "binomialProbabilities"),
  function(x) {
    length(x@grp1.prob)
  }
)


#' @export new.binomial
new.binomial <- function(x1, m1, x2, m2) {
  new("binomialInitialData",
    grp1.count = x1,
    grp2.count = x2,
    size = study.size(
      grp1.size = m1,
      grp2.size = m2
    )
  )
}


beta.sample <- function(nsim, event.count, group.size) {
  # sample from a beta distribution (a posterior using
  # Jeffreys' posterior)
  #
  # Args:
  #   nsim        : number of samples
  #   event.count : count of events observed
  #   group.size  : in a group of given size

  rbeta(nsim, event.count + 0.5, group.size - event.count + 0.5)
}


setMethod(
  "samplePosterior",
  signature(earlyStudy = "binomialInitialData", nsim = "numeric"),
  function(earlyStudy, nsim) {
    # Samples the event probabilities in two groups
    #
    # Args:
    #   x    : data from the initial study
    #   nsim : the number of simulations to do
    #
    #
    # Returns:
    #   A binomialProbabilities object giving the event probabilities
    #   in each group.

    grp1.prob <- beta.sample(
      nsim, earlyStudy@grp1.count,
      earlyStudy@size@grp1
    )
    grp2.prob <- beta.sample(
      nsim, earlyStudy@grp2.count,
      earlyStudy@size@grp2
    )
    new("binomialProbabilities",
      grp1.prob = grp1.prob, grp2.prob = grp2.prob
    )
  }
)


setMethod(
  "treatmentEffect",
  signature(posteriorSample = "binomialProbabilities"),
  function(posteriorSample) {
    posteriorSample@grp1.prob - posteriorSample@grp2.prob
  }
)

setClass("binomialStudy",
  contains = "twoArm",
  representation(
    direction = "integer",
    testFunction = "function"
  ),
  validity = function(object) {
    dirn <- object@direction
    if (length(dirn) != 1L) {
      "direction must be a scalar"
    } else if (!(dirn %in% c(-1L, 1L))) {
      "direction must be in {-1, 1}"
    } else {
      TRUE
    }
  }
)


##' make a binomial study
##' @export new.binomialStudy
new.binomialStudy <- function(size,
                              endpoint = c("cure", "mortality"),
                              method = c(
                                "normal", "chisq", "wald",
                                "bl", "ac"
                              ),
                              significance = as.numeric(NA),
                              hurdle = as.numeric(NA),
                              margin = as.numeric(NA)) {
  ## first find our endpoint
  endpoint <- match.arg(endpoint)
  direction <- switch(endpoint,
    cure = 1L,
    mortality = -1L,
    stop("unknown 'endpoint'")
  )
  ## and what we're doing
  method <- match.arg(method)
  ## here's the test function
  testFunction <- switch(method,
    normal = .gaussianTest,
    chisq = .chiTest,
    wald = .ni.wald,
    bl = .ni.bl,
    ac = .ni.ac,
    stop("unknown test method")
  )
  ## check that the arguments are sane
  if (method %in% c("normal", "chisq")) {
    if (!missing(margin)) {
      stop("you have specified 'margin' for an equivalence study")
    }
    assert(is.numeric(hurdle) ||
      (is.positive.scalar(significance) &&
        (significance < 1)))
  } else {
    if (!missing(hurdle)) {
      stop("you have specified 'hurdle' for a noninferiority study")
    }
    if (missing(margin)) {
      stop("you need to specify 'margin' for a noninferiority study")
    }
    if (missing(significance)) {
      stop("you need to specific 'significance' for a noninferiority study")
    }
    hurdle <- margin
    assert(is.numeric(hurdle))
    assert(is.positive.scalar(significance) && (significance < 1))
  }
  new("binomialStudy",
    size = size,
    significance = significance,
    hurdle = hurdle,
    direction = direction,
    testFunction = testFunction
  )
}


setClass(
  "binomialLater",
  representation(
    grp1.count = "numeric",
    grp2.count = "numeric",
    study.defn = "binomialStudy"
  )
)

setMethod(
  "sampleLater",
  signature(
    posteriorSample = "binomialProbabilities",
    laterStudy = "twoArm"
  ),
  function(posteriorSample, laterStudy) {
    n1 <- laterStudy@size@grp1
    n2 <- laterStudy@size@grp2
    nsim <- length(posteriorSample)
    grp1.count <- rbinom(nsim, n1, posteriorSample@grp1.prob)
    grp2.count <- rbinom(nsim, n2, posteriorSample@grp2.prob)
    new("binomialLater",
      grp1.count = grp1.count,
      grp2.count = grp2.count,
      study.defn = laterStudy
    )
  }
)


.chiTest <- function(event.1, size.1, event.2, size.2, alpha, hurdle,
                     flip.direction) {
  d <- ifelse(flip.direction * event.1 > flip.direction * event.2,
    1, 0
  )
  p.1 <- event.1 / size.1
  p.2 <- event.2 / size.2

  check.significance(
    alpha,
    {
      chi2 <- qchisq(alpha, 1, lower.tail = FALSE)
      n <- size.1 + size.2
      stat <- (n * (event.1 * (size.2 - event.2) -
        (size.1 - event.1) * event.2)^2) / (event.1 + event.2) / (n - event.1 - event.2) / size.1 / size.2 * d
      stat > chi2
    }
  ) & check.hurdle(flip.direction * (p.1 - p.2), hurdle)
}

.gaussianTest <- function(event.1, size.1, event.2, size.2, alpha, hurdle,
                          flip.direction) {
  p.1 <- event.1 / size.1
  p.2 <- event.2 / size.2
  effect <- flip.direction * (p.1 - p.2)

  check.significance(
    alpha,
    {
      se <- sqrt(p.1 * (1 - p.1) / size.1 + p.2 * (1 - p.2) / size.2)
      effect > se * qnorm(alpha / 2, lower.tail = FALSE)
    }
  ) & check.hurdle(effect, hurdle)
}

## backend code for non-inferiority tests
#' @name internal
.ni.test <- function(effect, se, alpha, margin, direction) {
  scale <- qnorm(alpha / 2, lower.tail = FALSE)
  lower.limit <- effect - scale * se
  upper.limit <- effect + scale * se
  if (direction > 0) {
    ## cure: we want p.1 > p.2
    lower.limit > margin
  } else {
    ## mortality: we want p.1 < p.2
    upper.limit < margin
  }
}

# Wald noninferiority test
#' @name internal
.ni.wald <- function(event.1, size.1, event.2, size.2, ...) {
  p.1 <- event.1 / size.1
  p.2 <- event.2 / size.2
  se <- sqrt(p.1 * (1 - p.1) / size.1 + p.2 * (1 - p.2) / size.2)
  if (any(se == 0)) {
    stop("Wald test: se estimate is zero, use a different CI method")
  }
  .ni.test(p.1 - p.2, se, ...)
}

# Agresti & Caffo noninferiority test
#' @name internal
.ni.ac <- function(event.1, size.1, event.2, size.2, ...) {
  p.1 <- (event.1 + 1) / (size.1 + 2)
  p.2 <- (event.2 + 1) / (size.2 + 2)
  se <- sqrt(p.1 * (1 - p.1) / size.1 + p.2 * (1 - p.2) / size.2)
  .ni.test(p.1 - p.2, se, ...)
}

# Brown & Lis Jeffreys method
#' @name internal
.ni.bl <- function(event.1, size.1, event.2, size.2, ...) {
  p.1 <- (event.1 + 0.5) / (size.1 + 1)
  p.2 <- (event.2 + 0.5) / (size.2 + 1)
  se <- sqrt(p.1 * (1 - p.1) / size.1 + p.2 * (1 - p.2) / size.2)
  .ni.test(p.1 - p.2, se, ...)
}

setMethod(
  "testLater",
  signature(laterSample = "binomialLater"),
  function(laterSample) {
    study.def <- laterSample@study.defn
    alpha <- study.def@significance
    sizing <- study.def@size
    study.def@testFunction(laterSample@grp1.count,
      sizing@grp1,
      laterSample@grp2.count,
      sizing@grp2,
      alpha,
      study.def@hurdle,
      study.def@direction)
  }
)
