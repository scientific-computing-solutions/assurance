

#' Internal classes and methods for assurance package
#' 
#' Internal assurance functions.
#' 
#' These internal \pkg{assurance} functions are not usually called directly by
#' the user.
#' 
#' @aliases assurance-internal.Rd gaussianCommonSDParameters-class
#' gaussianCommonSDParameters gaussianDifferentSDParameters-class
#' gaussianDifferentSDParameters length length-method testLater
#' @keywords internal
#' @name internal
NULL





#' Compute Assurance for Multiple Correlated Endpoints
#' 
#' These generic functions provide the framework for assurance computations for
#' trials with multiple endpoints.
#' 
#' These generic functions tie together the whole assurance calculation for
#' multiple endpoints.  By default they do Monte-Carlo simulation using
#' \code{samplePosterior} (which samples from the posterior distribution),
#' \code{sampleLater} (which, given a sample from the posterior, simulates the
#' later study), and \code{sampleAndTestLater}.  \code{sampleAndTestLater} is
#' just the composition of \code{samplePosterior} and \code{testLater}.
#' 
#' \code{sampleAndTestLater} exists for the purposes of shortcut calculations
#' in which, given a sample from the posterior, the probability of later stage
#' success can be directly determined.
#' 
#' @param earlyStudy Data defining the initial study.
#' @param laterStudy Information about the later study.
#' @param nsim Number of Monte-Carlo simulations (maybe unnecessary).
#' @param posteriorSample A sample from our posterior distribution.
#' @param laterSample Simulations of the later trial.
#' @return \code{assurance.multiple.endpoints} returns a list of items:
#' NoAdjustment.Each returns (uncorrected) assurance values for each endpoint
#' separately. NoAdjustment.AtLeastOne returns the (uncorrected) assurance that
#' at least one endpoint will be statistically significant. NoAdjustmentAll
#' returns the (uncorrected) assurance that all endpoints will be statistically
#' significant. Currently no multiplicity adjustments are implemented.
#' 
#' \code{samplePosterior} returns an opaque object that contains samples from
#' the posterior.
#' 
#' \code{sampleLater} returns an opaque object that contains simulated results
#' for the later study.
#' 
#' \code{testLater} returns probability of success for the simulated later
#' study.
#' 
#' \code{sampleAndTestLater} returns probability of success for simulated later
#' study.
#' @author Paul Metcalfe, Mary Jenner, Daniel Dalevi
#' @name assurance.multiple.endpoints
#' @export assurance.multiple.endpoints
NULL


#' Compute assurance
#' 
#' Perform assurance computations for a clinical trial, comparing means,
#' proportions, or hazards.
#' 
#' Estimate the unconditional probability of a study achieving its desired goal
#' by averaging the power of the study over the prior distribution of the
#' treatment effect. This is termed `assurance' by O'Hagan et al (2005), but
#' also see Spiegelhalter et al (200X) and the references therein.
#' 
#' The information provided to the function is used to construct a prior
#' density. For Gaussian data, the marginal density of the variance is taken to
#' be scaled inverse \eqn{\chi^2}{chi^2} and the conditional density of the
#' mean is taken to be Gaussian (see pages 74 - 75 of Gelman et al).
#' 
#' For binomial data, the marginal proportions are taken to be beta
#' distributions resulting from using Jeffreys' prior. That is, \eqn{p_1 \sim
#' B(x_1 + 0.5, m_1 - x_1 + 0.5)}{p1 ~ Beta(x1 + 0.5, m1 - x1 + 0.5} and
#' similarly for \eqn{p_2}{p2}.
#' 
#' For time-to-event data, the log of the hazard ratio is taken to be Gaussian
#' with known variance \eqn{1/x_1 + 1/x_2}{1/x1 + 1/x2} (see, for example,
#' Carroll, 2003).  Treating the variance as fixed fails to account for some
#' uncertainty, but in practice the loss is small because it is usual for the
#' number of events to be prespecified as part of the study design.
#' 
#' @name assurance-package
#' @docType package
#' @author Harry Southworth, Paul Metcalfe
#' 
#' Maintainer: Paul Metcalfe <paul.metcalfe@@astrazeneca.com>
#' @references A. O'Hagan, J. W. Stevens and M. J. Campbell, Assurance in
#' clinical trial design, Pharmaceutical Statistics,4, 187 - 201, 2005
#' 
#' D. J. Spiegelhalter, K. R. Abrams and J. P. Myles, Bayesian Approaches to
#' Clinical Trials and Health-care Evaluation, Wiley, 2003
#' 
#' A. Gelman, J. B. Carlin, H. S. Stern and D. B. Rubin, Bayesian Data Analysis
#' (Second Edition), Chapman & Hall/CRC, 2004
#' 
#' K. J. Carroll, On the use and utility of the Weibull model in the analysis
#' of survival data, Controlled Clinical Trials, 24, 682 - 701, 2003
#' @keywords design
#' @import methods stats4
#' @importFrom Rcpp sourceCpp
#' @useDynLib assurance, .registration=TRUE
NULL





#' Compute Assurance
#' 
#' These generic functions provide the framework for assurance computations.
#' 
#' These generic functions tie together the whole assurance calculation.  By
#' default they do Monte-Carlo simulation using \code{samplePosterior} (which
#' samples from the posterior distribution), \code{sampleLater} (which, given a
#' sample from the posterior, simulates the later study), and
#' \code{sampleAndTestLater}.  \code{sampleAndTestLater} is just the
#' composition of \code{samplePosterior} and \code{testLater}.
#' 
#' \code{sampleAndTestLater} exists for the purposes of shortcut calculations
#' in which, given a sample from the posterior, the probability of later stage
#' success can be directly determined.
#' 
#' @aliases assurance samplePosterior sampleLater testLater sampleAndTestLater
#' @param earlyStudy Data defining the initial study.
#' @param laterStudy Information about the later study.
#' @param nsim Number of Monte-Carlo simulations (maybe unnecessary).
#' @param posteriorSample A sample from our posterior distribution.
#' @param laterSample Simulations of the later trial.
#' @return \code{assurance} returns a scalar; the calculated assurance.
#' 
#' \code{samplePosterior} returns an opaque object that contains samples from
#' the posterior.
#' 
#' \code{sampleLater} returns an opaque object that contains simulated results
#' for the later study.
#' 
#' \code{testLater} returns probability of success for the simulated later
#' study.
#' 
#' \code{sampleAndTestLater} returns probability of success for simulated later
#' study.
#' @author Paul Metcalfe
#' @examples
#' 
#' count.data <-
#'   new.overdispersedPoisson(grp1.count=27,
#'                            grp1.overdispersion=1.6,
#'                            grp1.time.at.risk=97,
#'                            grp2.count=56,
#'                            grp2.overdispersion=1.6,
#'                            grp2.time.at.risk=82)
#' 
#' later.study <-
#'   new.overdispersedPoissonStudy(duration=56./52,
#'                                 grp1.size=228,
#'                                 grp2.size=228,
#'                                 grp1.overdispersion=1.6,
#'                                 grp2.overdispersion=1.6,
#'                                 significance=0.04)
#' 
#' assurance(count.data, later.study, nsim=100000)
#' @export assurance
#' @name assurance
NULL





#' Binomial assurance tests
#' 
#' Creates an object representing results of initial study.
#' 
#' 
#' @aliases new.binomial samplePosterior,binomialInitialData,numeric-method
#' sampleLater,binomialProbabilities,twoArm-method
#' testLater,binomialLater-method length,binomialProbabilities-method
#' @param x1 The number of `successes' observed in treatment group 1 in the
#' existing study.
#' @param m1 The number of observations from which \code{x1} successes were
#' observed
#' @param x2 The number of `successes' observed in treatment group 2 in the
#' existing study.
#' @param m2 The number of observations from which \code{x2} successes were
#' observed
#' @author Paul Metcalfe
#' @examples
#' 
#' bin.data <- new.binomial(x1=18, m1=20, x2=15, m2=20)
#' @name bin
NULL
