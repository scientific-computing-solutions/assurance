% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assuranceTable.R
\name{assuranceTable}
\alias{assuranceTable}
\title{Tabulated assurance results}
\usage{
assuranceTable(earlyStudy, laterStudy, nsim, ...)
}
\arguments{
\item{earlyStudy}{Data defining the initial study}

\item{laterStudy}{Data defining the later study}

\item{nsim}{Number of simulations of later trial}

\item{...}{Miscellaneous other arguments}
}
\value{
Returns a data frame with

\item{low}{Lower interval bounds.}

\item{high}{Upper interval bounds.}

\item{interval.probability}{The probability of getting a treatment effect in
the given bound, based on the early study data.}

\item{prob.of.success}{The probability of getting a positive result in the
later study.  May have NaN value if no simulations fall in given interval.}
}
\description{
Tabulate assurance results for presentation.
}
\examples{

  initial.data <- new.gaussian(delta.mu=6, sd1=10, sd2=20, m1=50, m2=20)
  later.study <- new.twoArm(size=study.size(grp1.size=100,
                                            grp2.size=200),
			    significance=0.05)
  assuranceTable(initial.data, later.study, 100000)
}
\author{
Paul Metcalfe
}
