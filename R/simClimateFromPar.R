#################################

parManager <- function(parS, SWGparameterization, datInd) {
  UseMethod("parManager", parS)
}

SWGsim <- function(SWGpar, nTimes, randomTerm, obs=NULL) {
  UseMethod("SWGsim", SWGpar)
}

#################################

simClim = function(parS,               # vector of pars (will change in optim)
                   modelTag,
                   ppTypes=NULL,
                   datInd,
                   randomTerm=NULL,
                   obs=NULL){
  
  SWGmodel = strsplit(modelTag,'-')[[1]][[3]]
  SWGparameterization = strsplit(modelTag,'-')[[1]][[2]]
  
  timeStep = modelInfoList[[modelTag]]$timeStep
  if (is.null(timeStep)){
    timeStep = obs$timeStep
  }
  
  names(parS) = modelInfoList[[modelTag]]$parNam
  class(parS) = SWGmodel
  
  SWGpar = parManager(parS = parS, SWGparameterization = SWGparameterization, datInd=datInd[[aggNameShort[[timeStep]]]])
  
  class(SWGpar) = SWGmodel
  
  sim = SWGsim(SWGpar = SWGpar,
                  nTimes = datInd[[aggNameShort[[timeStep]]]]$nTimes,
                  randomTerm = randomTerm,
                  obs=obs)

  simVar = modelInfoList[[modelTag]]$simVar
  for (pp in ppTypes){
    sim = runPP(sim,obs[[simVar]],pp,parS,datInd) 
  }
  
  return(sim)  
  
}

#################################
