plot_attributes_vs_perturbed_OAT = function(clim,sim,attsEval,attPerturb,...){

  P = calcPerformanceAttributes(clim=clim,sim=sim,attSel=attsEval)

  for (att in names(P)){
    plotPerformanceOAT(P, sim, metric=att,plotType='base',attSel=attPerturb,...)
  }
  
}



