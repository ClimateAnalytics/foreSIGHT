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

modelInfoList[['Temp-annWD-wgenLM']] = list(simVar='Temp',
                                          timeStep = '1 day',
                                          simPriority=2,
                                          npars=3,
                                          parNam = c("cor0","mu.W","sigma.W",
                                                     "mu.D","sigma.D"),
                                         minBound=c(0.45,-10,0.1,-10,1),
                                         maxBound=c(0.95,40,20,40,20),
                                          # minBound=c(0.8,10,5,20,10), 
                                          # maxBound=c(0.8,10,5,20,10),
                                          WDcondition=TRUE,
                                          WDthresh=0)

modelInfoList[['Temp-seasWD-wgenLM']] = list(simVar='Temp',
                                            timeStep = '1 day',
                                            simPriority=2,
                                            npars=20,
                                            parNam = c("cor0.SON","cor0.DJF","cor0.MAM","cor0.JJA",
                                                       "muW.SON","muW.DJF","muW.MAM","muW.JJA",
                                                       "sigmaW.SON","sigmaW.DJF","sigmaW.MAM","sigmaW.JJA",
                                                       "muD.SON","muD.DJF","muD.MAM","muD.JJA",
                                                       "sigmaD.SON","sigmaD.DJF","sigmaD.MAM","sigmaD.JJA"
                                            ),
                                            minBound=c(0.45,0.45,0.45,0.45,
                                                       -10,-10,-10,-10,
                                                       0.1,0.1,0.1,0.1,
                                                       10,10,10,10,
                                                       1,1,1,1),
                                            maxBound=c(0.95,0.95,0.95,0.95,
                                                       40,40,40,40,
                                                       20,20,20,20,
                                                       40,40,40,40,
                                                       20,20,20,20),
                                            # minBound=c(0.8,10,5,20,10), 
                                            # maxBound=c(0.8,10,5,20,10),
                                            WDcondition=TRUE,
                                            WDthresh=0,
                                            dependentVar='P')

modelInfoList[['PET-seasWD-wgenLM']] = list(simVar='PET',
                                             timeStep = '1 day',
                                             simPriority=2,
                                             npars=20,
                                             parNam = c("cor0.SON","cor0.DJF","cor0.MAM","cor0.JJA",
                                                        "muW.SON","muW.DJF","muW.MAM","muW.JJA",
                                                        "sigmaW.SON","sigmaW.DJF","sigmaW.MAM","sigmaW.JJA",
                                                        "muD.SON","muD.DJF","muD.MAM","muD.JJA",
                                                        "sigmaD.SON","sigmaD.DJF","sigmaD.MAM","sigmaD.JJA"
                                             ),
                                             minBound=c(0.45,0.45,0.45,0.45,
                                                        0,0,0,0,
                                                        0.1,0.1,0.1,0.1,
                                                        0,0,0,0,
                                                        0.1,0.1,0.1,0.1),
                                             maxBound=c(0.95,0.95,0.95,0.95,
                                                        10,10,10,10,
                                                        2,2,2,2,
                                                        10,10,10,10,
                                                        2,2,2,2),
                                             WDcondition=TRUE,
                                             WDthresh=0,
                                             dependentVar='P',
                                            minVal=0)

modelInfoList[['PET-har-wgenLM']] = list(simVar='PET',
                                          timeStep = '1 day',
                                          simPriority=2,
                                          npars=7,
                                          parNam = c("cor0",
                                                     "mu.m","mu.amp","mu.ang",
                                                     "sigma.m","sigma.amp","sigma.ang"),
                                          minBound=c(0.1,
                                                     0,0,0,
                                                     0,0,0),
                                          # maxBound=c(0.95,
                                          #            10,5,6.28,
                                          #            2,1,6.28),
                                         maxBound=c(0.95,
                                                    10,5,15,
                                                    2,1,15),
                                         WDcondition=FALSE)


modelInfoList[['PET-harWD-wgenLM']] = list(simVar='PET',
                                            timeStep = '1 day',
                                            simPriority=2,
                                            npars=13,
                                           parNam = c("cor0",
                                                      "muW.m","muW.amp","muW.ang",
                                                      "sigmaW.m","sigmaW.amp","sigmaW.ang",
                                                      "muD.m","muD.amp","muD.ang",
                                                      "sigmaD.m","sigmaD.amp","sigmaD.ang"),
                                            minBound=c(0.1,
                                                       0,0,-3.14,
                                                       0,0,-3.14,
                                                       0,0,-3.14,
                                                       0,0,-3.14),
                                            maxBound=c(0.95,
                                                       10,5,3.14,
                                                       2,1,3.14,
                                                       10,5,3.14,
                                                       5,2,3.14),
                                            WDcondition=TRUE,
                                            WDthresh=0,
                                            dependentVar='P',
                                            minVal=0)


modelInfoList[['PET-seas-wgenLM']] = list(simVar='PET',
                                           timeStep = '1 day',
                                           simPriority=2,
                                           npars=12,
                                           parNam = c('cor0.SON','cor0.DJF','cor0.MAM','cor0.JJA',
                                                      'mu.SON','mu.DJF','mu.MAM','mu.JJA',
                                                      'sigma.SON','sigma.DJF','sigma.MAM','sigma.JJA'),
                                           minBound=c(0.45,0.45,0.45,0.45,
                                                      0,0,0,0,
                                                      0.1,0.1,0.1,0.1), 
                                           maxBound=c(0.9,0.9,0.9,0.9,
                                                      10,10,10,10,
                                                      2,2,2,2),
                                           WDcondition=FALSE,
                                          minVal=0)



# #################################

parManager.wgenLM = function(parS, SWGparameterization, datInd, auxInfo=NULL){

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

  }   else if (SWGparameterization=='annWD'){
    cor0 = rep(parS['cor0'],datInd$nTimes)
    # change=c(1,which(abs(diff(auxInfo$wdStatus))==1)+1)
    # cor0[change] = 0
    mu = rep(parS['mu.D'],datInd$nTimes)
    mu[auxInfo$wdStatus] = parS['mu.W']
    sigma = rep(parS['sigma.D'],datInd$nTimes)
    sigma[auxInfo$wdStatus] = parS['sigma.W']

  }   else if (SWGparameterization=='seasWD'){
    cor0 <- assignSeasPars(parS['cor0.SON'], parS['cor0.DJF'], parS['cor0.MAM'], parS['cor0.JJA'], datInd[["i.ss"]])
    change=c(1,which(abs(diff(auxInfo$wdStatus))==1)+1)
    cor0[change] = 0
    
    muW <- assignSeasPars(parS['muW.SON'], parS['muW.DJF'], parS['muW.MAM'], parS['muW.JJA'], datInd[["i.ss"]])
    sigmaW <- assignSeasPars(parS['sigmaW.SON'], parS['sigmaW.DJF'], parS['sigmaW.MAM'], parS['sigmaW.JJA'], datInd[["i.ss"]])
    muD <- assignSeasPars(parS['muD.SON'], parS['muD.DJF'], parS['muD.MAM'], parS['muD.JJA'], datInd[["i.ss"]])
    sigmaD <- assignSeasPars(parS['sigmaD.SON'], parS['sigmaD.DJF'], parS['sigmaD.MAM'], parS['sigmaD.JJA'], datInd[["i.ss"]])
    
    mu = muD
    mu[auxInfo$wdStatus] = muW[auxInfo$wdStatus] 
    
    sigma = sigmaD
    sigma[auxInfo$wdStatus] = sigmaW[auxInfo$wdStatus] 
    
  }   else if (SWGparameterization=='harWD'){
    cor0 = rep(parS['cor0'],datInd$nTimes)
    change=c(1,which(abs(diff(auxInfo$wdStatus))==1)+1)
    cor0[change] = 0
  
    muW = harmonicFunc(x=1:datInd$nTimes,mean=parS['muW.m'],amp=parS['muW.amp'],phase.ang=parS['muW.ang'],k=1,nperiod=365)
    sigmaW = harmonicFunc(x=1:datInd$nTimes,mean=parS['sigmaW.m'],amp=parS['sigmaW.amp'],phase.ang=parS['sigmaW.ang'],k=1,nperiod=365)
    
    muD = harmonicFunc(x=1:datInd$nTimes,mean=parS['muD.m'],amp=parS['muD.amp'],phase.ang=parS['muD.ang'],k=1,nperiod=365)
    sigmaD = harmonicFunc(x=1:datInd$nTimes,mean=parS['sigmaD.m'],amp=parS['sigmaD.amp'],phase.ang=parS['sigmaD.ang'],k=1,nperiod=365)
    
    mu = muD
    mu[auxInfo$wdStatus] = muW[auxInfo$wdStatus] 
    
    sigma = sigmaD
    sigma[auxInfo$wdStatus] = sigmaW[auxInfo$wdStatus] 
    
  }
  
  parTS = list(cor0=cor0,mu=mu,sigma=sigma)
  
  return(parTS)
  
}


#################################

SWGsim.wgenLM = function(SWGpar,
                            nTimes,
                            randomTerm,
                            auxInfo=NULL){
  
  if (!is.null(randomTerm$randomUnitNormalVector)){
    randomUnitNormalVector = randomTerm$randomUnitNormalVector
  } else {
    randomUnitNormalVector=stats::qnorm(randomTerm$randomVector)
  }

  epsilonT = randomUnitNormalVector*SWGpar$sigma * sqrt(1-SWGpar$cor0^2)
  X = latentX_calc_cpp(SWGpar$cor0, epsilonT, nTimes)
  X = X + SWGpar$mu

  sim = X
      
  minVal = auxInfo$modelInfo$minVal
  if(!is.null(minVal)){
    sim[sim<minVal]=minVal
  }
  
  return(sim)
  
}

#################################
