#' @include assurance.R
NULL

##' Calculates assurace for survival data, assumes time-event data is a function of arm (either 0 or 1)
##' @param total.events Number of events we est. will happen in later study
##' @param alpha Significance level (default 0.05)
##' @param repl Number of simulations used in assurance calculation
##' @param C_ARM Column name of the arm data which is represented by 0 or 1 for active or control
##' @param C_TIME Column name of the event times
##' @param C_EVENT Column name of the indicator for the event 
##' @return A function that calculates the assurance values for survival data 
##' @export surv.statistic
surv.statistic = function( total.events, alpha=0.05, repl=1e5, C_ARM="TRTNUM", C_TIME="D_DAYS", C_EVENT="DEATH" ) {
  first = TRUE
  return( function( x, i ) {
      if( first ){ check.arm( x, C_ARM ); first = FALSE }
      re.data = data.frame( x[[ C_ARM ]][i], x[[ C_EVENT ]][i], x[[ C_TIME ]][i])
      names( re.data ) = c( "arm", "event", "time" )
      mod = coxph( Surv( time, event ) ~ arm, data=re.data ) 
      hz = exp( coef( mod ) )
      x1 = sum( re.data[ re.data$arm==0, ]$event ) 
      x2 = sum( re.data[ re.data$arm==1, ]$event ) 
      st.data = new.survival( hz, x1 = x1 , x2 = x2 ) 
      later.st = new.twoArm( size=study.size( total.size=total.events ), significance=alpha )
      assurance( st.data, later.st, repl )
 }) 
}

##' Calculates assurace for gaussian data (difference in means)
##' @param n0 Subjects in control arm in planned study
##' @param n1 Subjects in active arm in planned study
##' @param alpha Significance level (default 0.05)
##' @param repl Number of simulations used in assurance calculation
##' @param C_ARM Column name of the arm data which is represented by 0 or 1 for active or control
##' @param C_EFFECT Column name for the effects (endpoint values of interest)
##' @return A function that calculates the assurance values for gaussian data
##' @export gauss.statistic
gauss.statistic = function( n0, n1, alpha=0.05, repl=1e5, C_ARM="arm", C_EFFECT="effect" ){
  first = TRUE
  return( function( x, i ) {
    if( first ){ check.arm( x, C_ARM ); first = FALSE }
    re.data = data.frame( x[[ C_ARM ]][i], x[[ C_EFFECT ]][i] )
    names( re.data ) = c( "arm", "effect" )
    my.i = tapply( re.data$effect, re.data$arm, mean )
    sd.i = tapply( re.data$effect, re.data$arm, sd )
    delta.my = my.i[1] - my.i[2]
    x0 = sum( re.data$arm==0 ) 
    x1 = sum( re.data$arm==1 ) 
    st.data = new.gaussian( delta.mu=delta.my, sd1=sd.i[1], sd2=sd.i[2], m1=x0, m2=x1 )
    later.st =  new.twoArm( size=study.size( grp1.size=n0, grp2.size=n1 ), significance=alpha )
    assurance( st.data, later.st, repl )
  })
}

##' Calculates assurace for binomial data (patient having a reponse or
##' not) comparing two arms
##' @param n Total number of subjects in new study
##' @param alpha Significance level (default 0.05)
##' @param repl Number of simulations used in assurance calculation
##' @param C_ARM Column name of the arm data which is represented by 0 or 1 for active or control
##' @param C_EFFECT Column name for the effects (endpoint values of interest)
##' @return A function that calculates the assurance values for gaussian data
##' @export binomial.statistic
binomial.statistic = function( n, alpha=0.05, repl=1e5, C_ARM="arm", C_EFFECT="response", hurdle=NULL ) {
  first = TRUE
  return( function( x, i ) {
    if( first ){ check.arm( x, C_ARM ); first = FALSE }
    re.data = data.frame( x[[ C_ARM ]][i], x[[ C_EFFECT ]][i] )
    names( re.data ) = c( "arm", "response" )
    x.0 = sum( re.data[ re.data$arm==0, ]$response ) 
    x.1 = sum( re.data[ re.data$arm==1, ]$response ) 
    n.0 = sum( re.data$arm==0 ) 
    n.1 = sum( re.data$arm==1 )
    initial.data = new.binomial( x1=x.0,  m1=n.0, x2=x.1, m2=n.1 )
    if( is.null(hurdle) ) {
      later.study  = new.binomialStudy( size=study.size( total.size=n ), significance=alpha )
      return( assurance( initial.data, later.study, repl ) )
    }
    else {
      later.study  = new.binomialStudy( size=study.size( total.size=n ), significance=alpha, hurdle=hurdle )
      return( assurance( initial.data, later.study, repl ) )
    }
  })
}

##' Function that will check the ARM data column for integrity
check.arm = function( data, C_ARM ){
  if( sum( names(data) == C_ARM ) < 1 ){
    stop( paste0( "No column named C_ARM = ", C_ARM ) )
  }
  if( class(data[[C_ARM]]) != "factor" || sum(as.integer(levels(data[[C_ARM]]))) != 1 ) {
    stop("Column ARM need to be a factor with levels (0,1)")
  }
}


##' Function for bootstrapping assurance values
##' @param data Data frame containing the data to be
##' bootstrapped. Will contain different columns depending on which
##' type of analysis we wish to analyse.
##' @param statistic Function that will perform the assurance value
##' calculation
##' @param R number of bootstrap replications (Default: 500)
##' @export assurance.boot
##' @importFrom boot boot
##' @importFrom ggplot2 qplot geom_vline theme annotate element_text
assurance.boot = function( data,
                           statistic = surv.statistic( total.events=500 ),
                           R=500 ){ 
  if( sum( is.na( data ) ) > 0  ) {
    warning( "Dataset contains missing values. If these are in effect columns bootstrap will not proceed!" )
  }
  
  message("RUNNING BOOTSTRAP, ", R, " REPLICATIONS, PLEASE WAIT!")
  
  x.boot = boot( data, statistic, R )
  ci = as.numeric( quantile( x.boot$t, c( 0.1, 0.9 ) )  )
  options( warn = -1 )
    p = qplot( x.boot$t, geom=c( "histogram" ), xlab="Assurance", 
               main=paste0( "Bootstrap ",  R, " replications"), xlim=c(0,1)  ) + 
      geom_vline( xintercept=( x.boot$t0 ), color="red", lwd=1 ) + 
      annotate( "text", x=0.5,y=-R/50, size=5, color="black",
                label=paste0( "Obs=", x.boot$t0, ", Quantiles(10,90)=[", ci[1], ", ", ci[2], "]" )   ) + 
      theme( plot.title = element_text(size = 15, face = "bold", colour = "black", vjust = 1))
  suppressMessages( print(p) )
  options( warn = 0 )
  return( x.boot )
}

##' Plotting Kaplan Meier plot for survival data, used in examples.
##' @importFrom survival Surv coxph survfit
##' @export
kaplan.meier.plot = function( data, C_ARM="TRTNUM", C_TIME="D_DAYS", C_EVENT="DEATH" ){
    mod = coxph( Surv( data[[ C_TIME ]], data[[ C_EVENT ]] ) ~ data[[ C_ARM ]]) 
    hz = exp( coef( mod ) )
    plot( survfit(  Surv( data[[ C_TIME ]], data[[ C_EVENT ]] ) ~ data[[ C_ARM ]]), 
          main=paste0( "Kaplan-Meier plot, HAZARD=",hz ), 
          xlab="Time", 
          ylab="Events")
}

##' Function that will simulate a gaussian trial, used in examples.
##' @param n0 number of subjects in arm 0
##' @param n1 number of subjects in arm 1 
##' @param my0 mean effect in control
##' @param my1 mean effect in active
##' @param sd0 sd in control
##' @param sd1 sd in active
##' @export simulate.gauss.trial
simulate.gauss.trial = function( n0, n1, my0, my1, sd0, sd1 ) {
  arm1 = rnorm( n0, my0, sd0 )
  arm2 = rnorm( n1, my1, sd1 )
  arm = as.factor( c( rep( 0, n0 ), rep( 1, n1 ) ) )
  id = 1:(n0+n1)
  effect = c( arm1, arm2 )
  data.frame( id, arm, effect )
}

##' Simulate a gaussian trial, used in examples.
##' @param n0 number of subjects in arm 0
##' @param n1 number of subjects in arm 1 
##' @param p0 event probability in control
##' @param p1 event probability in active
##' @export simulate.binomial.trial
simulate.binomial.trial = function( n0, n1, p0, p1, max0=1, max1=1 ) {
  arm0 = rbinom( n0, max0, p0 )
  arm1 = rbinom( n1, max1, p1 )
  arm = as.factor( c( rep( 0, n0 ), rep( 1, n1 ) ) )
  id = 1:(n0+n1)
  response = c( arm0, arm1 )
  data.frame( id, arm, response )
}
