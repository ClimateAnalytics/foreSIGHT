# # scaling
# 
# settings.R
# 
# load data camels
# 
# expSpace
# 
# sim scaling
# 
# system response
# 
# plot OAT, 2d
# 
# 
# # stoch
# 
# settings.R
# 
# load data camels
# 
# expSpace - OAT (Ptot, P99, DSD) (hold - nwet, seas, ann CV)
# 
# model settings
# 
# sim stoch P (not shuffle)
# sim stoch P shuffle
# sim scaling PET
# 
# combine
# 
# system response
# 
# plot OAT


rm(list=ls())

source('testing_2.0/settings.R')
#source('settings.R')

#catchment = 'A5030502' # Scott Creek
catchment = 'A5050517' # North Para River at Penrice (A5050517)
startYr = 1976
endYr = 2005

data = load_camels(catchment,startYr,endYr)

# create reference data 
clim_ref <- list(times = data$times,
                 P = data$P,
                 PET = data$PET)


# specify perturbed attributes
#attPerturb <- c("PET_day_all_avg", "P_day_all_tot","P_day_all_seasRatioMarMay") # note that updates to attribute manager allow customisation of months used for seasRatio 
attPerturb <- c("PET_day_all_avg", "P_day_all_tot","P_day_all_seasRatioMarAug") # note that updates to attribute manager allow customisation of months used for seasRatio 

# specify perturbation type and minimum-maximum ranges of the perturbed attributes
attPerturbType <- "regGrid"
#attPerturbSamp <- c(3,5,5)     
#attPerturbMin <- c(1, 0.7,0.7)
#attPerturbMax <- c(1.2, 1.3,1.3)
attPerturbSamp <- c(3,3,3)     
attPerturbMin <- c(1, 0.7,0.7)
attPerturbMax <- c(1.2, 1.3,1.3)
# create the exposure space using foreSIGHT
expSpace <- createExpSpace(attPerturb = attPerturb,
                           attPerturbSamp = attPerturbSamp,            # set to null
                           attPerturbMin = attPerturbMin,
                           attPerturbMax = attPerturbMax,
                           attPerturbType = attPerturbType,
                           attHold = NULL)

# generate simulations using simple scaling with seasonality 
sim <- generateScenarios(reference = clim_ref,             # input observed data
                         expSpace = expSpace,        # exposure space created by the user
                         controlFile = "scaling")    # using simple scaling with seasonality


source('testing_2.0/GR4J_funcs.R')
  
dates = as.Date(data$times)

Param = setup_cal_GR4J(dates = dates,data$P,data$PET,data$Qobs)

systemArgs = list(dates=dates,Param=Param)

sysOutSim = runSystemModel(sim=sim,systemModel = GR4J_wrapper,systemArgs = systemArgs,metrics = c('meanQ','P99','P25'))
  
sysOutClim = GR4J_wrapper(data = clim_ref, systemArgs = systemArgs,metrics = c('meanQ','P99','P25'))
  
plotPerformanceOAT(performance = sysOutSim, sim=sim, metric = 'meanQ')
plotPerformanceOAT(performance = sysOutSim, sim=sim, metric = 'P99')
plotPerformanceOAT(performance = sysOutSim, sim=sim, metric = 'P25')

plotPerformanceSpace(performance = sysOutSim, sim=sim, metric = 'meanQ',attX = 'P_day_all_tot',attY='P_day_all_seasRatioMarMay')


