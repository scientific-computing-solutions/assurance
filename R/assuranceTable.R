#' @include assurance.R
NULL


##' Find good break points
findBreaks <- function(x, breaks = "Sturges", nclass = NULL) {
  use.br <- !missing(breaks)
  if (use.br) {
    if (!missing(nclass)) {
      warning("'nclass' not used when 'breaks' is specified")
    }
  } else if (!is.null(nclass) && length(nclass) == 1) {
    breaks <- nclass
  }
  use.br <- use.br && (nB <- length(breaks)) > 1
  if (use.br) {
    breaks <- sort(breaks)
  } else {
    if (is.character(breaks)) {
      breaks <- match.arg(
        tolower(breaks),
        c(
          "sturges", "fd",
          "freedman-diaconis", "scott"
        )
      )
      breaks <- switch(breaks,
        sturges = nclass.Sturges(x),
        `freedman-diaconis` = ,
        fd = nclass.FD(x),
        scott = nclass.scott(x),
        stop("unknown 'breaks' algorithm")
      )
    } else if (is.function(breaks)) {
      breaks <- breaks(x)
    }
    if (!is.numeric(breaks) || !is.finite(breaks) || breaks < 1) {
      stop("invalid number of 'breaks'")
    }
    breaks <- pretty(range(x), n = breaks, min.n = 1)
    nB <- length(breaks)
    if (nB <= 1) {
      stop("findBreaks: pretty() error, breaks=", format(breaks))
    }
  }
  if (any(diff(breaks) <= 0)) {
    stop("'breaks' are not strictly increasing")
  }
  breaks
}



#' Tabulated assurance results
#'
#' Tabulate assurance results for presentation.
#'
#'
#' @param earlyStudy Data defining the initial study
#' @param laterStudy Data defining the later study
#' @param nsim Number of simulations of later trial
#' @param ... Miscellaneous other arguments
#' @return Returns a data frame with
#'
#' \item{low}{Lower interval bounds.}
#'
#' \item{high}{Upper interval bounds.}
#'
#' \item{interval.probability}{The probability of getting a treatment effect in
#' the given bound, based on the early study data.}
#'
#' \item{prob.of.success}{The probability of getting a positive result in the
#' later study.  May have NaN value if no simulations fall in given interval.}
#' @author Paul Metcalfe
#' @examples
#'
#' initial.data <- new.gaussian(delta.mu = 6, sd1 = 10, sd2 = 20, m1 = 50, m2 = 20)
#' later.study <- new.twoArm(
#'   size = study.size(
#'     grp1.size = 100,
#'     grp2.size = 200
#'   ),
#'   significance = 0.05
#' )
#' assuranceTable(initial.data, later.study, 100000)
#' @export
assuranceTable <- function(earlyStudy, laterStudy, nsim, ...) {
  ## first simulate the treatment effect
  posterior <- samplePosterior(earlyStudy, nsim)
  assert(length(posterior) == nsim, "not enough samples")
  ## and now get the treatment effect
  effect <- treatmentEffect(posterior)
  breaks <- c(-Inf, findBreaks(effect, ...), Inf)
  num.breaks <- length(breaks)
  num.intervals <- num.breaks - 1
  intervals <- rbind(
    low = breaks[1:num.intervals],
    high = breaks[2:num.breaks]
  )
  probs <-
    apply(
      intervals, 2,
      function(col, effect, success) {
        low.bnd <- col[[1]]
        high.bnd <- col[[2]]
        flag <- (effect > low.bnd) & (effect <= high.bnd)
        prob.in.range <- mean(flag)
        power <- mean(success[flag])
        c(prob.in.range, power)
      }, effect, sampleAndTestLater(posterior, laterStudy)
    )
  rownames(probs) <- c("interval.probability", "prob.of.success")
  as.data.frame(t(rbind(intervals, probs)))
}
