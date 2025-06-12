######################

systemModel_calcAtts <- function(data,          # data.frame with columns: year, month, day, *var1*, *var2* etc.
                                 systemArgs,    # list containing the arguments of simulateSystem
                                 metrics) {     # names of performance metrics (with units of the metrics)
  
  if (!is.null(systemArgs$vSel)){
    if (systemArgs$cSel=='mean'){
      data[[systemArgs$vSel]] = apply(data[[systemArgs$vSel]],1,mean)
    } else {
      data[[systemArgs$vSel]] = data[[systemArgs$vSel]][,systemArgs$cSel]
    }
  }
  
  systemPerformanceSim = calculateAttributes(climateData=data,attSel=metrics)
  
  systemPerformance = (systemPerformanceSim/systemArgs$attBase-1)*100
  
  return(systemPerformance)
}

#################################################################

#' @export
calcPerformanceAttributes = function(clim,sim,attSel,vSel=NULL,cSel=1){

  if (!is.null(vSel)){
    if (cSel=='mean'){
      clim[[vSel]] = apply(clim[[vSel]],1,mean)
    } else {
      clim[[vSel]] = clim[[vSel]][,cSel]
    }
  }

  attBase = calculateAttributes(clim,attSel)
  
  systemPerf <- runSystemModel(sim = sim,                     # simulation; the perturbed time series
                               systemModel = systemModel_calcAtts,      # the system model function
                               systemArgs = list(attBase=attBase,
                                                 vSel=vSel,
                                                 cSel=cSel),        # argument to the system model function
                               metrics = attSel)              # selected performance metrics
  
  return(systemPerf)
  
}

##################################################

plotPerformanceAttributes = function(clim,sim,attPerturb,attEval){
  
  P = calcPerformanceAttributes(clim=tank_obs,sim=sim,attSel=attEval)

  attSel='P_day_all_tot_m'
  par(mfrow=c(4,4),mar=c(4,4,2,1))
  for (att in names(P)){
    plotPerformanceOAT(P, sim, metric=att,col='black',use_ggplot = F,attSel=attPerturb)
  }

}

##################################################