summarize_observations = function(obs=NULL, offsets=NULL, family="binomial", tweak=0.05) {

  if (family=="binomial") {
    # for binomial, prob=0,1, becomes infinite so minor fix for hyper parameter prior approximation
    # tweak = 0.05 # tail truncation prob
    obs[ obs==1 ] = 1 - tweak
    obs[ obs==0 ] = tweak

    obs = inla.link.logit( obs )
    if (!is.null(offsets)) {
      obs = obs - log( offsets )
    }

  } else if (family=="gaussian") {
    # do nothing
  } else {
    message( "Use of other families requires more parameterizations here .. ") 
    stop()
  }

  ll = which(is.finite(obs))

  mq = quantile( obs[ll], probs=c(0.025, 0.5, 0.975) )

  MS = c( 
    mean=mean(obs[ll]), 
    sd=sd(obs[ll]), 
    min=min(obs[ll]), 
    max=max(obs[ll]),  
    lb=mq[1], 
    median=mq[2], 
    ub=mq[3]  
  )  # on data /user scale not internal link

  return (MS)
}
