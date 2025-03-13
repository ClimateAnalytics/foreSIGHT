#rm(list=ls())

devtools::load_all()

######################

systemModel_calcAtts <- function(data,          # data.frame with columns: year, month, day, *var1*, *var2* etc.
                                 systemArgs,    # list containing the arguments of simulateSystem
                                 metrics) {     # names of performance metrics (with units of the metrics)
  
  #  climTmp = clim
  #  climTmp$P = data$P
  
  print('calculating P metrics')
  
  attSel = unlist(strsplit(metrics, split='Del_', fixed=TRUE))[seq(2,(2*length(metrics)),2)]
  
  systemPerformanceSim <- calculateAttributes(climateData=data,attSel=attSel,attCalcInfo=systemArgs$attCalcInfo)
  
  systemPerformanceBase = systemArgs$attBase
  
  systemPerformance = ( systemPerformanceSim - systemPerformanceBase) / systemPerformanceBase * 100
  
  systemPerformance[attSel=='P_ann_tot_cor'] = systemPerformanceSim[attSel=='P_ann_tot_cor']
  
  return(systemPerformance)
}

#################################################################

calcPerformance = function(clim,sim,attSel,baseTar=NULL){
  
  # tarName = paste0('Target',baseTar)
  # climBase = clim
  # nReps = sum(grepl('Rep',names(sim)))
  # attBase = 0
  # for (r in 1:nReps){
  #   repName = paste0('Rep',r)
  #   if(all(sim$controlFile == 'scaling')){
  #     climBase$P = sim[[repName]][[tarName]]$P
  #   } else {
  #     climBase$P = sim[[repName]][[tarName]]$P$sim
  #   }
  #   attBase = attBase + calculateAttributes(climBase,attSel)
  # }
  # attBase = attBase/nReps
  
  attCalcInfo = calculateAttributes(clim,attSel,return_attCalcInfo = T)
  attBase = calculateAttributes(clim,attSel,attCalcInfo=attCalcInfo)
  
  systemPerf = list()
  
  metrics = paste0('Del_',attSel)
  
  #  metrics = attSel
  
  systemPerf <- runSystemModel(sim = sim,                     # simulation; the perturbed time series
                               systemModel = systemModel_calcAtts,      # the system model function
                               systemArgs = list(attBase=attBase,
                                                 attCalcInfo=attCalcInfo),        # argument to the system model function
                               metrics = metrics)              # selected performance metrics
  
  return(systemPerf)
  
}

##################################################

plotPerformanceOATdm = function(attVal,
                                metricList,
                                metricName,
                                ylim=NULL,xlim=NULL,
                                y_normalized=T,
                                attPerturb,
                                colLine='black',colShade='grey',lwd=1.5,add=F,
                                base='rep'){
  
  yMat = metricVal[[metric]]
  
  if (is.null(dim(yMat))){
    yMat = matrix(yMat,nrow=1)
  }
  
  nTar = dim(yMat)[1]
  
  x = attVal
  # #  if (attPerturb != 'P_ann_tot_cor'){
  # if (!attPerturb%in%c('P_ann_tot_cor','P_ann_tot_corLong')){
  #   x = (x-1)*100
  # }
  
  if (y_normalized){
    b = which(x==0)
    #b = which(abs(x)==min(abs(x)))
    if (base=='med'){
      yBase = median(yMat[b,])
      yMat = (yMat/yBase-1)*100
    } else if (base=='rep'){
      yBase = yMat[b,]
      for (t in 1:nTar){
        yMat[t,] = (yMat[t,]/yBase-1)*100
      }      
    }
  }
  
  if (dim(yMat)[2]>1){
    med = apply(yMat,1,median)
    p5 = apply(yMat,1,quantile,0.05)
    p95 = apply(yMat,1,quantile,0.95)
  } else {
    med = yMat
  }
  
  if (!add){
    plot(x,med,type='l',xlim=xlim,ylim=ylim,xaxs='i',xlab='',ylab='',col=colLine)
    abline(h=0,lty=2,lwd=0.5,col='darkgrey')
    abline(v=0,lty=2,lwd=0.5,col='darkgrey')
  }
  
  if (dim(yMat)[2]>1){
    polygon(c(rev(x), x), c(rev(p95), p5), col = colShade, border = NA)
  }
  lines(x,med,type='l',col=colLine,lwd=lwd) 
  box()
  
  points(x,med,col=colLine,lwd=lwd) 
  
}

##################################################

# #attPerturb<-c("P_ann_tot_m","P_ann_seasRatio")
# attPerturb<-c("P_day_all_tot_m","P_day_all_seasRatio")
# #attPerturb<-c("P_year_all_avg","P_day_all_seasRatio")
# attPerturbType = "regGrid"
# attPerturbSamp = c(9, 9)
# attPerturbMin = c(0.8, 0.8)
# attPerturbMax = c(1.2, 1.2)
# expSpace <- createExpSpace(attPerturb = attPerturb,
#                            attPerturbSamp = attPerturbSamp,
#                            attPerturbMin = attPerturbMin,
#                            attPerturbMax = attPerturbMax,
#                            attPerturbType = attPerturbType)
# data(tankDat)
# seasScaling <- generateScenarios(reference = tank_obs,
#                                  expSpace = expSpace,
#                                  controlFile = "scaling")
# # load(file='testing_2.0/sim_seas_LV.RData')
# # sim = sim_seas

##################################################

# Example 3: Stochastic simulation using foreSIGHT default settings
#----------------------------------------------------------------------
attPerturb <- c("P_day_all_tot_m")
attHold <- c("P_day_all_nWet_m","P_day_Feb_tot_m", "P_day_SON_dyWet_m", "P_day_JJA_avgWSD_m", "P_day_MAM_tot_m",
             "P_day_DJF_avgDSD_m")
attPerturbType = "regGrid"
attPerturbSamp = c(5)
attPerturbMin = c(0.8)
attPerturbMax = c(1.2)

expSpace <- createExpSpace(attPerturb = attPerturb,
                           attPerturbSamp = attPerturbSamp,
                           attPerturbMin = attPerturbMin,
                           attPerturbMax = attPerturbMax,
                           attPerturbType = attPerturbType,
                           attHold = attHold)
# load example data available in foreSIGHT
data(tankDat)
# perform stochastic simulation
simStochastic <- generateScenarios(reference = tank_obs,
                                   expSpace = expSpace,
                                   seedID = 1,
                                   numReplicates = 5)


plotScenarios(simStochastic)

pause

##################################################

sim = seasScaling

attSel = c(attPerturb,attHold)

P = calcPerformance(clim=clim,sim=sim,attSel=attSel)

for (metric in names(P)){
  plotPerformanceSpace(P, sim, metric=metric)
  plotPerformanceOAT(P, sim, metric=metric)
  
}

plotPerformanceOATdm(attVal = P, sim, metric=metric)
