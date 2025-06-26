#' @export
evaluate_system_metrics = function(sim,clim,systemModel,systemArgs,metrics,obs_metrics=NULL,varNames=NULL){
  
  expSpace = sim$expSpace
  
  #########
  # determine target number that represent baseline unperturbed climate 
  attSel = colnames(sim$expSpace$targetMat)
  varType <- vapply(attSel, FUN = get.attribute.varType, FUN.VALUE = character(1), USE.NAMES = FALSE)
  targetType <- vapply(varType, FUN = get.target.type, FUN.VALUE = character(1), USE.NAMES = FALSE)
  
  baseVal = rep(NA,length(attSel))
  baseVal[targetType=='diff'] = 0
  baseVal[targetType=='frac'] = 1
  
  targets = expSpace$targetMat
  b = which(apply(targets==baseVal,1,FUN =all))
  
  #########
  # strip other targets from sim
  simBase = sim
  numReps = length(which(grepl('Rep',names(sim))))
  numTars = length(names(sim[[1]]))
  for (r in 1:numReps){
    repName = paste0('Rep',r)
    simBase[[repName]] = NULL
    simBase[[repName]][['Target1']] = sim[[repName]][[paste0('Target',b)]] 
  }
  
  #########
  # performance for baseline unperturbed climate  
  systemPerf_base <- runSystemModel(sim = simBase,                     # simulation; the perturbed time series
                                    systemModel = systemModel,      # the system model function
                                    systemArgs = systemArgs,        # argument to the system model function
                                    metrics = metrics,
                                    varNames=varNames)              # selected performance metrics 
  
  # performance using observed climate 
  systemPerf_obsClim = systemModel(data = clim, systemArgs = systemArgs, metrics = metrics)
  
  #############
  
  return(list(systemPerf_base=systemPerf_base,
              systemPerf_obsClim=systemPerf_obsClim))
  
}
