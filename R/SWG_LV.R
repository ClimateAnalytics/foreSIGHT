#################################

modelInfoList[["P-ann-LV"]] = list(simVar="P",
                                   timeStep = "1 day",
                                   simPriority=1,
                                   npars=4,
                                   parNam=c("alpha", "sigma", "mu", "lambda"),
                                   minBound=c(0, 0.001, -15, 1),
                                   maxBound=c(0.999, 10, 0, 2))
                                   # minBound=c(0, 0.001, -15, 0.5),
                                   # maxBound=c(0.999, 10, 5, 4))

modelInfoList[["P-seas-LV"]] = list(simVar="P",
                                    timeStep = "1 day",
                                    simPriority=1,
                                    npars=16,
                                    parNam=c("alpha.SON","alpha.DJF","alpha.MAM","alpha.JJA",
                                             "sigma.SON","sigma.DJF","sigma.MAM","sigma.JJA",
                                             "mu.SON","mu.DJF","mu.MAM","mu.JJA",
                                             "lambda.SON","lambda.DJF","lambda.MAM","lambda.JJA"),
                                    minBound=c(0,0,0,0,
                                               0.001,0.001,0.001,0.001,
                                               -15,-15,-15,-15,
                                               1,1,1,1),
                                    maxBound=c(0.999,0.999,0.999,0.999,
                                               10,10,10,10,
                                               0,0,0,0,
                                               4,4,4,4))


modelInfoList[['P-har-LV']] = list(simVar='P',
                                   timeStep = '1 day',
                                   simPriority=1,
                                   npars=12,
                                   parNam = c('alpha.m','alpha.amp','alpha.ang',
                                                'sigma.m','sigma.amp','sigma.ang',
                                                'mu.m','mu.amp','mu.ang',
                                                'lambda.m','lambda.amp','lambda.ang'),
                                   minBound=c(0, 0, 0,
                                              0.001, 0, 0,
                                              -15, 0, 0,
                                              1, 0, 0),
                                   maxBound=c(0.999, 0, 0,
                                              10, 5, 6.28,
                                              0, 8, 6.28,
                                              2, 0, 0))

# #################################

parManager.LV = function(parS, SWGparameterization, datInd){
  
  if (SWGparameterization=='ann'){
    alpha <- rep(parS['alpha'],datInd$nTimes)
    sigma <- rep(parS['sigma'],datInd$nTimes)
    mu <- rep(parS['mu'],datInd$nTimes)
    lambda <- rep(parS['lambda'],datInd$nTimes)
  } else if (SWGparameterization=='seas'){
    alpha <- assignSeasPars(parS['alpha.SON'], parS['alpha.DJF'], parS['alpha.MAM'], parS['alpha.JJA'], datInd[["i.ss"]])
    sigma <- assignSeasPars(parS['sigma.SON'], parS['sigma.DJF'], parS['sigma.MAM'], parS['sigma.JJA'], datInd[["i.ss"]])
    mu <- assignSeasPars(parS['mu.SON'], parS['mu.DJF'], parS['mu.MAM'], parS['mu.JJA'], datInd[["i.ss"]])
    lambda <- assignSeasPars(parS['lambda.SON'], parS['lambda.DJF'], parS['lambda.MAM'], parS['lambda.JJA'], datInd[["i.ss"]])
  } else if (SWGparameterization=='har'){
    alpha = harmonicFunc(x=seq(1:datInd$nTimes),mean=parS['alpha.m'],amp=parS['alpha.amp'],phase.ang = parS['alpha.ang'],k=1,nperiod=365)
    sigma = harmonicFunc(x=seq(1:datInd$nTimes),mean=parS['sigma.m'],amp=parS['sigma.amp'],phase.ang = parS['sigma.ang'],k=1,nperiod=365)
    mu = harmonicFunc(x=seq(1:datInd$nTimes),mean=parS['mu.m'],amp=parS['mu.amp'],phase.ang = parS['mu.ang'],k=1,nperiod=365)
    lambda = harmonicFunc(x=seq(1:datInd$nTimes),mean=parS['lambda.m'],amp=parS['lambda.amp'],phase.ang = parS['lambda.ang'],k=1,nperiod=365)
  }
  
  parTS = list(alpha=alpha,sigma=sigma,mu=mu,lambda=lambda)
  
  return(parTS)
  
}

#################################

SWGsim.LV = function(SWGpar,
                        nTimes,
                        randomTerm,
                        obs=NULL){
  
  if (length(SWGpar[['alpha']])!=nTimes){
    stop('length alpha != nTimes')
  }
  if (length(SWGpar[['sigma']])!=nTimes){
    stop('length sigma != nTimes')
  }
  if (length(SWGpar[['mu']])!=nTimes){
    stop('length mu != nTimes')
  }
  if (length(SWGpar[['lambda']])!=nTimes){
    stop('length lambda != nTimes')
  }
  
  if (!is.null(randomTerm$randomUnitNormalVector)){
    randomUnitNormalVector = randomTerm$randomUnitNormalVector
  } else {
    randomUnitNormalVector=stats::qnorm(randomTerm$randomVector)
  }
  
  if (length(randomUnitNormalVector)!=nTimes){
    stop('length randomUnitNormalVector != nTimes')
  }
  
  # Calculate latent variable - latentX
  epsilonT = randomUnitNormalVector*SWGpar$sigma
  X = latentX_calc_cpp(SWGpar$alpha, epsilonT, nTimes)
  X = X + SWGpar$mu
  
  rain = rep(0,nTimes)
  rain[X>0] = X[X>0] ^ SWGpar$lambda[X>0]

  return(rain)
  
}

#################################
