check_sim = function(sim,tol=0.05){

  nRep = length(which(grepl('Rep',names(sim))))
  nTarget = length(sim[['Rep1']])
  
  for (iRep in 1:nRep) {
    for (iTarg in 1:nTarget) {
      expTarg <- expSpace
      expTarg$targetMat <- expSpace$targetMat[iTarg, ]
      if(!is.null(expSpace$attRot)) {
        expTarg$attRot <- expSpace$attRot[iTarg]
      }
      if (!is.null(sim[[iRep]][[iTarg]]$attSim)){ # check if stochastic simulation performed
        varNames = names(sim[[iRep]][[iTarg]])
        varNames = varNames[!varNames%in%c('attSim','targetSim','parS','score')]
        
        for (var in varNames){
          
          if (any(sim[[iRep]][[iTarg]][[var]]$onBounds)){
            cat(paste0('parameters for ', var,' stoch rep, ', iRep, ' for target ',iTarg, ' on bounds\n'))
          }
          
          targDiff = abs(sim[[iRep]][[iTarg]]$targetSim - expTarg$targetMat) 
          if (any(targDiff > tol)){
            cat(paste0('error in target atts for ', var,' stoch rep ', iRep, ' for target ',iTarg, ' greater than tol\n'))
          }
          
        }
        
      }
      
    }
  }
  
}