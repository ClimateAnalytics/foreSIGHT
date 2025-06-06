rm(list=ls())

devtools::load_all()

#load('testing_2.0/test_plot_diagnostics.RData')

###################
# seasonal scaling

# attPerturb<-c("P_day_all_tot_m","P_day_all_seasRatioMarMay")
# attPerturbType = "regGrid"
# attPerturbSamp = c(5, 5)
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


##################################################

#Example 3: Stochastic simulation using foreSIGHT default settings
#----------------------------------------------------------------------
attPerturb = c("P_day_all_tot_m")
attHold = c("P_day_all_nWet_m","P_day_Feb_tot_m", "P_day_SON_dyWet_m", "P_day_JJA_avgWSD_m", "P_day_MAM_tot_m",
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

save.image(file='testing_2.0/test_plot_diagnostics.RData')


##################################################

# sim = seasScaling
sim = simStochastic

# attSel = c(attPerturb,attHold)

attOther = c('P_day_all_tot_m','P_day_SON_tot_m','P_day_DJF_tot_m','P_day_MAM_tot_m','P_day_JJA_tot_m',
             'P_day_all_P99','P_day_SON_P99','P_day_DJF_P99','P_day_MAM_P99','P_day_JJA_P99')

attSel = unique(c(attPerturb,attHold,attOther))
# attSel = unique(c(attPerturb,attOther))

P = calcPerformanceAttributes(clim=tank_obs,sim=sim,attSel=attSel)

attSel='P_day_all_tot_m'
par(mfrow=c(4,4),mar=c(4,4,2,1))
for (att in names(P)){
  # plotPerformanceSpace(P, sim, metric=metric)
  plotPerformanceOAT(P, sim, metric=att,col='black',use_ggplot = F,attSel=attSel)
}

attSel='P_day_all_seasRatioMarMay'
par(mfrow=c(4,4),mar=c(4,4,2,1))
for (att in names(P)){
  plotPerformanceOAT(P, sim, metric=att,col='black',use_ggplot = F,attSel=attSel)
}

names(P) = paste0('.',names(P))

attSel='P_day_all_tot_m'
#par(mfrow=c(4,4),mar=c(4,4,2,1))
for (metric in names(P)){
  tar=sim$expSpace$targetMat
  att = strsplit(metric,'[.]')[[1]][2]
  if (att%in%names(tar)){
    tar[metric] = (tar[att]-1)*100
  }
  plotPerformanceSpace(P, sim, metric=metric,climData = tar,axesPercentLabel='percentage.change')
}

  