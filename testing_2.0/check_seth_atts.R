rm(list=ls())

devtools::load_all()

#load('C:/Users/a1065639/Work/foreSIGHT/testing_2.0/seth.out_V2.RData')
#load('C:/Users/a1065639/Work/foreSIGHT/testing_2.0/seth.out_V2a.RData')
#load('C:/Users/a1065639/Work/foreSIGHT/testing_2.0/seth.out_V2b.RData')
#load('C:/Users/a1065639/Work/foreSIGHT/testing_2.0/seth.out_V2c.RData')
load('C:/Users/a1065639/Work/foreSIGHT/testing_2.0/seth.out_V2d.RData')

plotScenarios(sim)

attEval = colnames(sim$expSpace$targetMat)

attEval = c(attEval,
            'P_day_DJF_P99','P_day_MAM_P99','P_day_JJA_P99','P_day_SON_P99')
attEval = unique(attEval)

attEval = sort(attEval)

#attEval = attEval[1:5]

Perf = calcPerformanceAttributes(clim=clim_ref,sim=sim,attSel=attEval)

par(mfrow=c(5,6),mar=c(3,3,1,1))
attPerturb = 'P_day_all_P99'
plotPerformanceAttributesOAT(sim=sim,clim = clim_ref,attPerturb = attPerturb,attEval = attEval,Perf=Perf)

pause

par(mfrow=c(5,6),mar=c(3,3,1,1))
attPerturb = 'P_day_all_seasRatioDecFeb'
plotPerformanceAttributesOAT(sim=sim,clim = clim_ref,attPerturb = attPerturb,attEval = attEval,Perf=Perf)

pause

Perf1 = Perf
names(Perf1) = paste0('att_',names(Perf1))
for (metric in names(Perf1)){
  Perf_ave = apply(Perf1[[metric]],1,mean)
  zAbsMax = max(abs(Perf_ave),10)
  zlim = c(-zAbsMax,zAbsMax)
  
  plotPerformanceSpace(performance = Perf1,sim=sim,metric=metric,colLim=zlim,
                       colMap = RColorBrewer::brewer.pal(11, "PiYG"),type = 'heat.plot')
  
  # plotPerformanceSpace(performance = Perf1,sim=sim,metric=metric,colLim=zlim,
  #                      colMap = RColorBrewer::brewer.pal(11, "PiYG"),type = 'filled.contour')
  

  # plotPerformanceSpace(performance = Perf1,sim=sim,metric=metric,colLim=zlim,
  #                      colMap = colorRampPalette(c("#8E0152","#C51B7D","#DE77AE","#F1B6DA",
  #                                                  "#FDE0EF","#F7F7F7","#E6F5D0","#B8E186",
  #                                                  "#7FBC41","#4D9221","#276419")),
  #                      type = 'filled.contour')

    
  
  
  # plotPerformanceSpace(performance = Perf1,sim=sim,metric=metric,colLim=zlim,type = 'filled.contour')
  
}
