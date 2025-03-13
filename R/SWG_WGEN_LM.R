# #################################
 
#' @include default_parameters.R


modelInfoList[['Temp-ann-wgenLM']] = list(simVar='Temp',
                                          timeStep = '1 day',
                                          simPriority=2,
                                          npars=3,
                                          parNam = c("cor0","mu","sigma"),
                                          minBound=c(0.45,-10,1), 
                                          maxBound=c(0.9,40,20),
                                          WDcondition=FALSE)

modelInfoList[['Temp-seas-wgenLM']] = list(simVar='Temp',
                                              timeStep = '1 day',
                                              simPriority=2,
                                              npars=12,
                                              parNam = c('cor0.SON','cor0.DJF','cor0.MAM','cor0.JJA',
                                                           'mu.SON','mu.DJF','mu.MAM','mu.JJA',
                                                           'sigma.SON','sigma.DJF','sigma.MAM','sigma.JJA'),
                                              minBound=c(0.45,0.45,0.45,0.45,
                                                         -10,-10,-10,-10,
                                                         1,1,1,1,1), 
                                              maxBound=c(0.9,0.9,0.9,0.9,
                                                         40,40,40,40,
                                                         20,20,20,20,20),
                                              WDcondition=FALSE)

modelInfoList[['Temp-har-wgenLM']] = list(simVar='Temp',
                                       timeStep = '1 day',
                                       simPriority=2,
                                       npars=7,
                                       parNam = c("cor0",
                                                    "mu.m","mu.amp","mu.ang",
                                                    "sigma.m","sigma.amp","sigma.ang"),
                                       minBound=c(0.45,7.0,1.0,-0.05,0.9,0.1,-1.6),
                                       maxBound=c(0.9,28.0,9.0,0.81,4.9,1.4,3.15),
                                       WDcondition=FALSE)

# #################################

parManager.wgenLM = function(parS, SWGparameterization, datInd){

  if (SWGparameterization=='ann'){
    
    cor0 = rep(parS['cor0'],datInd$nTimes)
    mu = rep(parS['mu'],datInd$nTimes)
    sigma = rep(parS['sigma'],datInd$nTimes)
    
  } else if (SWGparameterization=='har'){

    cor0 = rep(parS['cor0'],datInd$nTimes)
    mu = harmonicFunc(x=1:datInd$nTimes,mean=parS['mu.m'],amp=parS['mu.amp'],phase.ang=parS['mu.ang'],k=1,nperiod=365)
    sigma = harmonicFunc(x=1:datInd$nTimes,mean=parS['sigma.m'],amp=parS['sigma.amp'],phase.ang=parS['sigma.ang'],k=1,nperiod=365)
 
  }   else if (SWGparameterization=='seas'){
    cor0 <- assignSeasPars(parS['cor0.SON'], parS['cor0.DJF'], parS['cor0.MAM'], parS['cor0.JJA'], datInd[["i.ss"]])
    mu <- assignSeasPars(parS['mu.SON'], parS['mu.DJF'], parS['mu.MAM'], parS['mu.JJA'], datInd[["i.ss"]])
    sigma <- assignSeasPars(parS['sigma.SON'], parS['sigma.DJF'], parS['sigma.MAM'], parS['sigma.JJA'], datInd[["i.ss"]])
  }
  
  parTS = list(cor0=cor0,mu=mu,sigma=sigma)
  
  return(parTS)
  
}


#################################

SWGsim.wgenLM = function(SWGpar,
                            nTimes,
                            randomTerm,
                            obs=NULL){
  
  if (!is.null(randomTerm$randomUnitNormalVector)){
    randomUnitNormalVector = randomTerm$randomUnitNormalVector
  } else {
    randomUnitNormalVector=stats::qnorm(randomTerm$randomVector)
  }
  
  
  epsilonT = randomUnitNormalVector*SWGpar$sigma
  X = latentX_calc_cpp(SWGpar$cor0, epsilonT, nTimes)
  X = X + SWGpar$mu
  
  sim = X
       
  return(sim)
  
}

#################################
