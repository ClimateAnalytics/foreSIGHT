#################################

#' @include default_parameters.R

modelInfoList[["P-ann-BLRPM"]] = list(simVar="P",
                                      timeStep='1 hour',
                                      simPriority=1,
                                      npars=5,
                                      parNam=c("lambda","gamma","beta","eta","mux"),
                                      # minBound=c(0.01,0.05,0.2,1.,3),
                                      # maxBound=c(0.03,0.15,0.4,3.,5))
                                      minBound=c(0.001,0.01,0.05,0.1,0.1),
                                      maxBound=c(0.1,0.5,1,10,10))

#################################

parManager.BLRPM = function(parS, SWGparameterization, datInd,auxInfo=NULL){
  
  if (SWGparameterization=='ann'){
    lambda <- parS['lambda']
    gamma <- parS['gamma']
    beta <- parS['beta']
    eta <- parS['eta']
    mux <- parS['mux']
  }
  # } else if (SWGparameterization=='seas'){
  #   mu <- assignSeasPars(parS['mu.SON'], parS['mu.DJF'], parS['mu.MAM'], parS['mu.JJA'], datInd[["i.ss"]])
  #   sigma <- assignSeasPars(parS['sigma.SON'], parS['sigma.DJF'], parS['sigma.MAM'], parS['sigma.JJA'], datInd[["i.ss"]])
  #   phi <- assignSeasPars(parS['phi.SON'], parS['phi.DJF'], parS['phi.MAM'], parS['phi.JJA'], datInd[["i.ss"]])
  #   lambda <- assignSeasPars(parS['lambda.SON'], parS['lambda.DJF'], parS['lambda.MAM'], parS['lambda.JJA'], datInd[["i.ss"]])
  # }
  
  parTS = list(lambda=lambda,gamma=gamma,beta=beta,eta=eta,mux=mux)
  
  return(parTS)
  
}

#################################

SWGsim.BLRPM = function(SWGpar,
                        nTimes,
                        randomTerm,
                        auxInfo=NULL){
  
  set.seed(randomTerm$seed)
  
  P = BLRPM.sim(lambda=SWGpar$lambda,
                gamma=SWGpar$gamma,
                beta=SWGpar$beta,
                eta=SWGpar$eta,
                mux=SWGpar$mux,
                t.sim=nTimes)$RR
  
  return(P)
  
}

#################################
