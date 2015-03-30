#' @include assurance.R
NULL

##' Assurance for Gaussian endpoints
##'
##' Assurance can be calculated for Gaussian endpoints, either
##' assuming common std. dev. between arms or different
##' std. dev. between arms.
##'
##' The prior distributions used are the improper
##' translation-invariant prior for treatment effect and the improper
##' scale invariant prior for the nuisance parameters of standard
##' deviation.
##' 
##' @aliases length,gaussianCommonSDParameters-method
##' length,gaussianDifferentSDParameters-method
##' samplePosterior,gaussianCommonSDInitialData,numeric-method
##' samplePosterior,gaussianDifferenceSDInitialData,numeric-method
##' treatmentEffect,gaussianCommonSDParameters-method
##' treatmentEffect,gaussianDifferentSDParameters-method
##' sampleLater,gaussianCommonSDParameters,twoArm-method
##' sampleLater,gaussianDifferentSDParameters,twoArm-method
##' testLater,gaussianCommonSDLater-method
##' testLater,gaussianDifferentSDLater-method
##' @name gaussian
##' @examples
##' initial.data <- new.gaussian(delta.mu=6,
##'                              sd1=10, sd2=20, m1=50, m2=20)
##' later.study <- new.twoArm(size=study.size(grp1.size=100,
##'                          grp2.size=200),
##'                       significance=0.05)
##' assurance(initial.data, later.study, 100000)
NULL

setClass("gaussianCommonSDInitialData",
         representation(delta.mu="numeric",
                        sigma2="numeric",
                        size="twoWayStudySize"),
         validity=function(object) {
           if (length(object@delta.mu) != 1) {
             "delta.mu must be a scalar"
           } else if (!is.positive.scalar(object@sigma2)) {
             "sigma2 must be a positive scalar"
           } else {
             TRUE
           }
         })


setClass("gaussianDifferentSDInitialData",
         representation(delta.mu="numeric",
                        grp1.sigma2="numeric",
                        grp2.sigma2="numeric",
                        size="twoWayStudySize"),
         validity=function(object) {
           if (length(object@delta.mu) != 1) {
             "delta.mu must be a scalar"
           } else if (!is.positive.scalar(object@grp1.sigma2)) {
             "grp1.sigma2 must be a positive scalar"
           } else if (!is.positive.scalar(object@grp2.sigma2)) {
             "grp2.sigma2 must be a positive scalar"
           } else {
             TRUE
           }
         })


setClass( "gaussianMultipleEndpointsInitialData",
          representation( mu1 = "numeric",
                          mu2 = "numeric",
                          ss.mat1 = "matrix",
                          ss.mat2 = "matrix",
                          size="twoWayStudySize" ),
          validity=function(object) {
            if( length( object@mu1 ) < 1) {
              "mu1 must have at least one element"
            }
            else {
              TRUE
            }
          })



setClass("gaussianMultipleEndpointsSDParameters",
         representation( mu.control = "list",
                         mu.trt     = "list" ) )




setClass("gaussianCommonSDParameters",
         representation(delta.mu="numeric",
                        sigma2="numeric"),
         validity=function(object) {
           if (length(object@delta.mu) != length(object@sigma2)) {
             return("delta.mu and sigma2 must be the same length")
           }
           if (any(object@sigma2 <= 0)) {
             return("sigma2 must be non-negative")
           }
           TRUE
         })


setMethod("length",
          signature(x="gaussianCommonSDParameters"),
          function(x) {
            length(x@delta.mu)
          })



setClass("gaussianDifferentSDParameters",
         representation(delta.mu="numeric",
                        grp1.sigma2="numeric",
                        grp2.sigma2="numeric"),
         validity=function(object) {
           len1 <- length(object@delta.mu)
           len2 <- length(object@grp1.sigma2)
           len3 <- length(object@grp2.sigma2)
           if (!((len1==len2) && (len2 == len3))) {
             return("delta.mu, grp1.sigma2, and grp2.sigma2 must be the same length")
           }
           if (any(object@grp1.sigma2 <= 0)) {
             return("grp1.sigma2 must be nonnegative")
           }
           if (any(object@grp2.sigma2 <= 0)) {
             return("grp2.sigma2 must be nonnegative")
           }
           TRUE
         })


setMethod("length",
          signature(x="gaussianDifferentSDParameters"),
          function(x) {
            length(x@delta.mu)
          })


##' @export new.gaussian
new.gaussian <- function(delta.mu, m1, m2, sd=NULL, sd1=NULL, sd2=NULL ) {
  stopifnot( m1 > 0 && m2 > 0 )
  if (is.null(sd1) && is.null(sd2) && !is.null(sd)) {
    new("gaussianCommonSDInitialData",
        delta.mu=delta.mu,
        sigma2=sd*sd,
        size=study.size(grp1.size=m1,
          grp2.size=m2))

  }
  else if (!is.null(sd1) && !is.null(sd2) && is.null(sd)) {
    new("gaussianDifferentSDInitialData",
        delta.mu=delta.mu,
        grp1.sigma2=sd1*sd1,
        grp2.sigma2=sd2*sd2,
        size=study.size(grp1.size=m1,
          grp2.size=m2))
  } else {
    stop("bad combination of parameters in new.gaussian")
  }
}


##' @export new.mult.gaussian
new.mult.gaussian <- function( ss.mat1, ss.mat2, mu1, mu2, m1, m2)
{
  stopifnot( m1 > 0 && m2 > 0 )
  if( dim( ss.mat1 )[1] != dim( ss.mat2)[1] ){
    stop( "#rows/cols of cov-matrices must agree")
  }
  if( nrow( ss.mat1 ) != length( mu1 ) ) {
    stop( "#rows/cols of cov-matrices must agree with size of mu vector")
  }
  new( "gaussianMultipleEndpointsInitialData",
       mu1 = mu1,
       mu2 = mu2,
       ss.mat1 = ss.mat1,
       ss.mat2 = ss.mat2,
       size = study.size( grp1.size = m1,
                          grp2.size = m2 ))
}

# herewith a scaled inverse chi^2 that matches Gelman p574, p580.

rsichisq <- function(nsim, df, scale) {
  scale * df / rchisq(nsim, df)
}


setMethod("samplePosterior",
          signature(earlyStudy="gaussianCommonSDInitialData", nsim="numeric"),
          function(earlyStudy, nsim) {
            # sample from the posterior for Gaussian RVS with common SD
            # degrees of freedom is grp1.size + grp2.size - 2
            sizes <- earlyStudy@size
            dof <- sizes@grp1 + sizes@grp2 - 2
            # the posterior for sigma2 is scaled inverse chi^2 with
            # given degrees of freedom
            sigma2 <- rsichisq(nsim, dof, earlyStudy@sigma2)
            # the variance of the difference of two means
            # is as given
            mean.variance <- sigma2 * (1/sizes@grp1 +
                                       1/sizes@grp2)
            # so the sampled mean differences are
            mean.sample <- rnorm(nsim,
                                 mean=earlyStudy@delta.mu,
                                 sd=sqrt(mean.variance))
            new("gaussianCommonSDParameters",
                delta.mu=mean.sample,
                sigma2=sigma2)
          })


setMethod( "samplePosterior",
           signature( earlyStudy="gaussianMultipleEndpointsInitialData",
                      nsim="numeric" ),
           function( earlyStudy, nsim ) {
             sizes <- earlyStudy@size
             n = length( earlyStudy@mu1 )
             
             trt =     postmvn( ybar=earlyStudy@mu2, S=earlyStudy@ss.mat2, sample.n = sizes@grp2, nsim=nsim )
             control = postmvn( ybar=earlyStudy@mu1, S=earlyStudy@ss.mat1, sample.n = sizes@grp1, nsim=nsim )
             
             new( "gaussianMultipleEndpointsSDParameters",
                  mu.control = control,
                  mu.trt     = trt )
           })


setMethod("treatmentEffect",
          signature(posteriorSample="gaussianCommonSDParameters"),
          function(posteriorSample) {
            posteriorSample@delta.mu
          })

setMethod("samplePosterior",
          signature(earlyStudy="gaussianDifferentSDInitialData",
                    nsim="numeric"),
          function(earlyStudy, nsim) {
            sizes <- earlyStudy@size
            grp1.sigma2 <- rsichisq(nsim, sizes@grp1 - 1,
                                    earlyStudy@grp1.sigma2)
            grp2.sigma2 <- rsichisq(nsim, sizes@grp2 - 1,
                                    earlyStudy@grp2.sigma2)
            mean.variance <- (grp1.sigma2 / sizes@grp1 +
                              grp2.sigma2 / sizes@grp2)
            mean.sample <- rnorm(nsim,
                                 mean=earlyStudy@delta.mu,
                                 sd=sqrt(mean.variance))
            new("gaussianDifferentSDParameters",
                delta.mu=mean.sample,
                grp1.sigma2=grp1.sigma2,
                grp2.sigma2=grp2.sigma2)
          })

setMethod("treatmentEffect",
          signature(posteriorSample="gaussianDifferentSDParameters"),
          function(posteriorSample) {
            posteriorSample@delta.mu
          })

sample.sigma2 <- function(nsim, scale, dof) {
   scale * rchisq(nsim, dof) / dof
}


#jeffery's prior for posterior
#' @importFrom mvtnorm rmvnorm
#' @importFrom MCMCpack riwish
postmvn=function( ybar, S, sample.n, nsim=100000)
{
  n = length(ybar)
  
  Lambdan=S
  nun=sample.n-1
  
  SIG=lapply(1:nsim,
             function(x, df, ss) riwish( v=df, S=ss ),
             df = nun, ss = Lambdan)
  
  MEAN=t( sapply(SIG, function(s)rmvnorm(1,ybar,s/sample.n)) )
  
  COV <- matrix(unlist(SIG), nrow = nsim, byrow = TRUE)
  
  x=list(MEAN,COV)
  names(x)=c("MEAN","COV")
  x
}



setClass("gaussianCommonSDLater",
         representation(delta.mu.sample="numeric",
                        variance.sample="numeric",
                        study.defn="twoArm"))


setMethod("sampleLater",
          signature(posteriorSample="gaussianCommonSDParameters",
                    laterStudy="twoArm"),
          function(posteriorSample, laterStudy) {
            sizing <- laterStudy@size
            n1 <- sizing@grp1
            n2 <- sizing@grp2
            nsim <- length(posteriorSample)
            delta.mu.sample <- rnorm(nsim,
                                     posteriorSample@delta.mu,
                                     sqrt(posteriorSample@sigma2 * (1/n1 + 1/n2)))
            dof <- n1 + n2 - 2
            variance.sample <- sample.sigma2(nsim,
                                             posteriorSample@sigma2, dof)
            new("gaussianCommonSDLater",
                delta.mu.sample=delta.mu.sample,
                variance.sample=variance.sample,
                study.defn=laterStudy)
          })


setClass("gaussianDifferentSDLater",
         representation(delta.mu.sample="numeric",
                        grp1.sigma2.sample="numeric",
                        grp2.sigma2.sample="numeric",
                        study.defn="twoArm"))


setMethod("sampleLater",
          signature(posteriorSample="gaussianDifferentSDParameters",
                    laterStudy="twoArm"),
          function(posteriorSample, laterStudy) {
            sizing <- laterStudy@size
            n1 <- sizing@grp1
            n2 <- sizing@grp2
            nsim <- length(posteriorSample)
            delta.mu.sample <- rnorm(nsim,
                                     posteriorSample@delta.mu,
                                     sqrt(posteriorSample@grp1.sigma2 / n1 +
                                          posteriorSample@grp2.sigma2 / n2))
            grp1.sigma2.sample <- sample.sigma2(nsim,
                                                posteriorSample@grp1.sigma2,
                                                n1 - 1)
            grp2.sigma2.sample <- sample.sigma2(nsim,
                                                posteriorSample@grp2.sigma2,
                                                n2 - 1)
            new("gaussianDifferentSDLater",
                delta.mu.sample=delta.mu.sample,
                grp1.sigma2.sample=grp1.sigma2.sample,
                grp2.sigma2.sample=grp2.sigma2.sample,
                study.defn=laterStudy)
          })

setClass( "gaussianMultipleEndpointsSDLater",
         representation( mu.control="matrix",
                         mu.trt = "matrix", 
                         cov.control = "matrix",
                         cov.trt = "matrix",
                         test.type = "character",
                         study.defn="twoArm" ) )


# No need for this... Just propagate the results
setMethod( "sampleLater",
           signature( posteriorSample="gaussianMultipleEndpointsSDParameters",
                      laterStudy="twoArm" ),
          function( posteriorSample, laterStudy ) {
            new("gaussianMultipleEndpointsSDLater",
                mu.control  = posteriorSample@mu.control$MEAN,
                mu.trt      = posteriorSample@mu.trt$MEAN,
                cov.control = posteriorSample@mu.control$COV,
                cov.trt     = posteriorSample@mu.trt$COV,
                study.defn  = laterStudy )
          })


setMethod( "testLater",
          signature(laterSample="gaussianCommonSDLater"),
          function(laterSample) {
            alpha <- laterSample@study.defn@significance
            sizing <- laterSample@study.defn@size
            check.significance(
              alpha,
              {
                n1 <- sizing@grp1
                n2 <- sizing@grp2
                dof <- n1 + n2 - 2
                t <- laterSample@delta.mu.sample /
                  sqrt(laterSample@variance.sample * (1/n1 + 1/n2))
                crit <- qt(0.5 * alpha, lower.tail=FALSE, df=dof)
                t > crit
              }) & check.hurdle(laterSample@delta.mu.sample,
                                laterSample@study.defn@hurdle)
          })


setMethod( "testLater",
          signature(laterSample="gaussianDifferentSDLater"),
          function(laterSample) {
            alpha <- laterSample@study.defn@significance
            sizing <- laterSample@study.defn@size

            check.significance(
              alpha,
              {
                n1 <- sizing@grp1
                n2 <- sizing@grp2
                delta.mu.sample <- laterSample@delta.mu.sample
                grp1.sigma2.sample <- laterSample@grp1.sigma2.sample
                grp2.sigma2.sample <- laterSample@grp2.sigma2.sample

                trm1 <- grp1.sigma2.sample / n1
                trm2 <- grp2.sigma2.sample / n2
                t <- delta.mu.sample / sqrt(trm1 + trm2)
                dof.use <- (trm1 + trm2)**2 / (trm1**2 / (n1 - 1) +
                                               trm2**2 / (n2 - 1))
                crit <- qt(0.5 * alpha, lower.tail=FALSE, df=dof.use)
                t > crit
              }) & check.hurdle(laterSample@delta.mu.sample,
                                laterSample@study.defn@hurdle)
          })


uncorrected.significant = function( laterSample )
{
  sizing <- laterSample@study.defn@size
  alpha <- laterSample@study.defn@significance
  N = nrow( laterSample@mu.control )
  m = ncol( laterSample@mu.control )
  
  delta <- laterSample@mu.trt - laterSample@mu.control
  hurdle <- laterSample@study.defn@hurdle
  
  significant = matrix( 0, N, m )
  
  for( i in 1:N ) 
  {
    
    trt.mean = laterSample@mu.trt[ i, ]
    control.mean = laterSample@mu.control[ i, ]  
    
    trt.cov =     matrix( laterSample@cov.trt[i,], m, m )
    control.cov = matrix( laterSample@cov.control[i,], m, m )
    
    cov = ( laterSample@study.defn@m1*trt.cov + laterSample@study.defn@m2*control.cov ) /
      ( laterSample@study.defn@m1 + laterSample@study.defn@m2 )
    cov.diag = sqrt( diag( diag( cov ) ) )
    cov.inv = solve( cov.diag )
    test.corr=cov.inv %*% cov %*% cov.inv
    
    
    denom = sqrt( 1/sizing@grp1 + 1/sizing@grp2) * sqrt( diag( cov ) )
    test.mean = ( trt.mean - control.mean )/denom
    test = rmvnorm( n = 1, mean=test.mean, sigma=test.corr )
    pv = 1 - pnorm( test )
    significant[ i, ] = ( pv <= alpha/2 )
    
    if (!any(is.na(hurdle)))
    {
      significant[ i, ] = significant[i, ] & (delta[i,] > hurdle)
    }
  }
  
  
  
  each = apply( significant, 2, mean )
  atleastone = sum( apply( significant, 1, sum ) >= 1 )/N
  all = sum( apply( significant, 1, sum ) == m )/N
  res = list( each, atleastone, all  ) 
  names( res ) = c( "NoAdjustment.Each", "NoAdjustment.AtLeastOne", "NoAdjustment.All" )
  
  
  
  return( res )
}

setMethod( "testLater",
          signature( laterSample="gaussianMultipleEndpointsSDLater" ),
          function( laterSample ) {
            alpha <- laterSample@study.defn@significance
            sizing <- laterSample@study.defn@size

            if (laterSample@study.defn@test.me == "Uncorrected.Significant") { 
                check.significance( alpha, 
                                    { 
                                      uncorrected.significant( laterSample ) 
                                      })  
            }
    })

