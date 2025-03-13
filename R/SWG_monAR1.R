#################################

#' @include default_parameters.R

modelInfoList[["P-ann-monAR1"]] = list(simVar="P",
                                       timeStep='1 month',
                                       simPriority=1,
                                       npars=4,
                                       parNam=c("mu","sigma","phi","lambda"),
                                       # minBound=c(9,7.,0.2,1.2),
                                       # maxBound=c(11,8.,0.3,1.6))
minBound=c(-1e2,0.01,-0.7,0.5),
maxBound=c(5e2,5e2,0.9,2))

#parS_ann = c(10,7.5,0.25,1.4)


modelInfoList[["P-seas-monAR1"]] = list(simVar="P",
                                        timeStep='1 month',
                                        simPriority=1,
                                       npars=16,
                                       parNam=c("mu.SON","mu.DJF","mu.MAM","mu.JJA",
                                                "sigma.SON","sigma.DJF","sigma.MAM","sigma.JJA",
                                                "phi.SON","phi.DJF","phi.MAM","phi.JJA",
                                                "lambda.SON","lambda.DJF","lambda.MAM","lambda.JJA"),
                                       minBound=rep(c(-1e2,0.01,-0.5,0.5),each=4),
                                       maxBound=rep(c(5e2,1e3,0.9,2),each=4))


modelInfoList[["P-har-monAR1"]] = list(simVar="P",
                                       timeStep='1 month',
                                       simPriority=1,
                                       npars=12,
                                       parNam=c("mu.m","mu.amp","mu.ang",
                                                "sigma.m","sigma.amp","sigma.ang",
                                                "phi.m","phi.amp","phi.ang",
                                                "lambda.m","lambda.amp","lambda.ang"),
                                       minBound=c(-1e2, 0, 0,
                                                  0.01, 0, 0,
                                                  -0.5, 0, 0,
                                                  0.5, 0, 0),
                                       maxBound=c(2e2, 2e2, 6.28,
                                                  2e2, 2e2, 6.28,
                                                  0.9, 1, 6.28,
                                                  2, 1, 6.28))

# #################################

parManager.monAR1 = function(parS, SWGparameterization, datInd){

  if (SWGparameterization=='ann'){
    mu <- rep(parS['mu'],datInd$nTimes)
    sigma <- rep(parS['sigma'],datInd$nTimes)
    phi <- rep(parS['phi'],datInd$nTimes)
    lambda <- rep(parS['lambda'],datInd$nTimes)
  } else if (SWGparameterization=='seas'){
    mu <- assignSeasPars(parS['mu.SON'], parS['mu.DJF'], parS['mu.MAM'], parS['mu.JJA'], datInd[["i.ss"]])
    sigma <- assignSeasPars(parS['sigma.SON'], parS['sigma.DJF'], parS['sigma.MAM'], parS['sigma.JJA'], datInd[["i.ss"]])
    phi <- assignSeasPars(parS['phi.SON'], parS['phi.DJF'], parS['phi.MAM'], parS['phi.JJA'], datInd[["i.ss"]])
    lambda <- assignSeasPars(parS['lambda.SON'], parS['lambda.DJF'], parS['lambda.MAM'], parS['lambda.JJA'], datInd[["i.ss"]])
  } else if (SWGparameterization=='har'){
    mu = harmonicFunc(x=seq(1:datInd$nTimes),mean=parS['mu.m'],amp=parS['mu.amp'],phase.ang = parS['mu.ang'],k=1,nperiod=12)
    sigma = harmonicFunc(x=seq(1:datInd$nTimes),mean=parS['sigma.m'],amp=parS['sigma.amp'],phase.ang = parS['sigma.ang'],k=1,nperiod=12)
    phi = harmonicFunc(x=seq(1:datInd$nTimes),mean=parS['phi.m'],amp=parS['phi.amp'],phase.ang = parS['phi.ang'],k=1,nperiod=12)
    lambda = harmonicFunc(x=seq(1:datInd$nTimes),mean=parS['lambda.m'],amp=parS['lambda.amp'],phase.ang = parS['lambda.ang'],k=1,nperiod=12)
    
  }

  sigma[sigma<0] = 0.
  phi[phi<-0.9]=-0.9
  phi[phi>0.9]=0.9
  
  parTS = list(mu=mu,sigma=sigma,phi=phi,lambda=lambda)

  # browser()
  
  return(parTS)

}

#################################

SWGsim.monAR1 = function(SWGpar,
                         nTimes,
                         randomTerm,
                         obs=NULL){

  if (length(SWGpar[['mu']])!=nTimes){
    stop('length mu != nTimes')
  }
  if (length(SWGpar[['sigma']])!=nTimes){
    stop('length sigma != nTimes')
  }
  if (length(SWGpar[['phi']])!=nTimes){
    stop('length phi != nTimes')
  }
  if (length(SWGpar[['lambda']])!=nTimes){
    stop('length lambda != nTimes')
  }
  if (length(randomTerm$randomVector)!=nTimes){
    stop('length randomTerm != nTimes')
  }

  sigma = SWGpar[['sigma']]*sqrt(1-SWGpar[['phi']]^2)
  epsilonT = qnorm(randomTerm$randomVector) * sigma
   # epsilonT = qnorm(randomTerm$randomVector) * SWGpar[['sigma']]
  X = latentX_calc_cpp(SWGpar[['phi']], epsilonT, nTimes)

  X = X + SWGpar[['mu']]
  
  P = rep(0,nTimes)
  P[X>0] = X[X>0] ^ SWGpar$lambda[X>0]
  
  if(is.na(sum(P))){browser()}
  
  return(P)

}

#################################
