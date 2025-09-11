#################################

rm(list=ls())

devtools::load_all()

clim = convert_climYMD_POSIXct(tank_obs)

#################################


modelInfoList[["P-ann-latentX"]] = list(simVar="P",
                                       timeStep = "1 day",
                                       simPriority=1,
                                       npars=4,
                                       parNam=c("alpha", "sigma", "mu", "lambda"),
                                       minBound=c(0, 0.001, -15, 0.5),
                                       maxBound=c(0.999, 10, 5, 4))

# #################################

parManager.latentX = function(parS, SWGparameterization, datInd, auxInfo=NULL){
  
  if (SWGparameterization=='ann'){
    parTS = assignAnnualParameters(parS=parS,datInd=datInd) 
  } 
  return(parTS)
}

#################################

SWGsim.latentX = function(SWGpar,
                         nTimes,
                         randomTerm,
                         auxInfo=NULL){
  
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

####################

modelSelection = list()
modelSelection$modelType = list()
modelSelection$modelType$P = "latentX"
modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "ann"
modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")

write(modelSelectionJSON, file = controlFile)

####################

attPerturb = c("P_day_all_tot_m")
attHold = c('P_day_all_P99','P_day_all_nWet','P_day_all_cor')

attPerturbType = "regGrid"
attPerturbSamp = c(1)
attPerturbMin = c(1)
attPerturbMax = c(1)

# create the exposure space
expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType,
                          attHold = attHold)


####################

sim = generateScenarios(reference = clim,
                        expSpace = expSpace,
                        controlFile = controlFile,
                        seedID=1)

plotScenarios(sim)

