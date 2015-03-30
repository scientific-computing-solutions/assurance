#' @include assurance.R
NULL

setClass("hazardInitialData",
         representation(log.hazard.ratio="numeric",
                        size="twoWayStudySize"))

setClass("hazardRatios",
         representation(log.hazard.ratio="numeric"))


setMethod("length",
          signature(x="hazardRatios"),
          function(x) {
            length(x@log.hazard.ratio)
          })


##' @export new.survival
new.survival <- function(hazard.ratio, x1, x2) {
  new("hazardInitialData",
      log.hazard.ratio=log(hazard.ratio),
      size=study.size(grp1.size=x1, grp2.size=x2))
}


setMethod("samplePosterior",
          signature(earlyStudy="hazardInitialData", nsim="numeric"),
          function(earlyStudy, nsim) {
            x <- earlyStudy
            stdev <- sqrt(1/x@size@grp1 + 1/x@size@grp2)
            log.hazard <- rnorm(nsim,
                                mean=x@log.hazard.ratio,
                                sd=stdev)
            new("hazardRatios", log.hazard.ratio=log.hazard)
          })


setMethod("treatmentEffect",
          signature(posteriorSample="hazardRatios"),
          function(posteriorSample) {
            exp(posteriorSample@log.hazard.ratio)
          })


setClass("hazardLater",
         representation(log.hazard="numeric",
                        study.defn="twoArm"))


setMethod("sampleLater",
          signature(posteriorSample="hazardRatios",
                    laterStudy="twoArm"),
          function(posteriorSample, laterStudy) {
            sizing <- laterStudy@size
            n1 <- sizing@grp1
            n2 <- sizing@grp2
            sd <- sqrt(1/n1 + 1/n2)
            nsim <- length(posteriorSample)
            log.hazard <- rnorm(nsim,
                                mean=posteriorSample@log.hazard.ratio,
                                sd=sd)
            new("hazardLater",
                log.hazard=log.hazard,
                study.defn=laterStudy)
          })


setMethod("testLater",
          signature(laterSample="hazardLater"),
          function(laterSample) {
            alpha <- laterSample@study.defn@significance
            sizing <- laterSample@study.defn@size
            check.significance(
              alpha,
              {
                n1 <- sizing@grp1
                n2 <- sizing@grp2
                sd <- sqrt(1/n1 + 1/n2)
                test.statistic <- laterSample@log.hazard / sd
                crit <- -qnorm(0.5 * laterSample@study.defn@significance,
                               lower.tail=FALSE)
                test.statistic < crit
              }) & check.hurdle(-laterSample@log.hazard,
                                -laterSample@study.defn@hurdle)
          })

# we are able to do the calculation explicitly in this case!

hazard.variance <- function(x) {
  1/x@grp1 + 1/x@grp2
}

setMethod("assurance",
          signature(earlyStudy="hazardInitialData",
                    laterStudy="twoArm",
                    nsim="missing"),
          function(earlyStudy, laterStudy, nsim) {
            # see the vignette for details of this calculation
            initial.sigma2 <- hazard.variance(earlyStudy@size)
            secondary.sigma2 <- hazard.variance(laterStudy@size)
            effect.variance <- initial.sigma2 + secondary.sigma2

            crit <- -qnorm(0.5 * laterStudy@significance, lower.tail=FALSE)
            lhr <- earlyStudy@log.hazard.ratio
            pnorm(min(crit * sqrt(secondary.sigma2),
                      laterStudy@hurdle,
                      na.rm=TRUE),
                  mean=lhr,
                  sd=sqrt(effect.variance))
          })
