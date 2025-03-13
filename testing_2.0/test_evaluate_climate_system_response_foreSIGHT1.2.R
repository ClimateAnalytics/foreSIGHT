rm(list=ls())

devtools::load_all('C:/Users/a1065639/Work/foreSIGHT1.2/foreSIGHT/')

data(tankDat)

#############################
# calc baseVal from expSpace 

evaluate_system_metrics = function(sim,clim,systemModel,systemArgs,metrics,obs_metrics=NULL){
  
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
                                    metrics = metrics)              # selected performance metrics 
  
  # performance using observed climate 
  systemPerf_obsClim = systemModel(data = tank_obs, systemArgs = systemArgs, metrics = metrics)
  
  #############
  
  return(list(systemPerf_base=systemPerf_base,
              systemPerf_obsClim=systemPerf_obsClim))
  
}

#############################################################

# specify perturbed attributes
attPerturb <- c("Temp_ann_avg_m", "P_ann_tot_m")

# specify perturbation type and minimum-maximum ranges of the perturbed attributes
attPerturbType <- "regGrid"
#attPerturbSamp <- c(9, 13)
#attPerturbSamp <- c(3, 3)
#attPerturbMin = c(-1, 0.80)
#attPerturbMax = c(1, 1.2)

attPerturbSamp <- c(1, 1)
attPerturbMin = c(0, 1)
attPerturbMax = c(0, 1)

# create the exposure space
expSpace <- createExpSpace(attPerturb = attPerturb, 
                           attPerturbSamp = attPerturbSamp, 
                           attPerturbMin = attPerturbMin,
                           attPerturbMax = attPerturbMax, 
                           attPerturbType = attPerturbType,
                           attHold = NULL)                    # no attributes held at historical levels

# generate perturbed time series using simple scaling
# simScaling <- generateScenarios(reference = tank_obs,             # input observed data
#                          expSpace = expSpace,        # exposure space created by the user
#                          controlFile = "scaling")    # using simple scaling

#############################################################

attPerturb <- c("Temp_ann_avg_m", "P_ann_tot_m")
attHold <- c("P_Feb_tot_m", "P_ann_nWet_m", "P_ann_R10_m", "P_SON_dyWet_m",
             "P_JJA_avgWSD_m", "P_MAM_tot_m",
             "P_DJF_avgDSD_m", "Temp_ann_rng_m")
attPerturbSamp <- c(1, 1)
attPerturbMin = c(0, 1)
attPerturbMax = c(0, 1)
# attPerturb <- c("P_day_all_tot_m")
# attHold <- c("P_day_Feb_tot_m", "P_day_all_nWet_m", "P_day_all_R10_m", "P_day_SON_dyWet_m", 
#              "P_day_JJA_avgWSD_m", "P_day_MAM_tot_m",
#              "P_day_DJF_avgDSD_m")
# attPerturbSamp <- c(1)
# attPerturbMin = c(1)
# attPerturbMax = c(1)
attPerturbType = "regGrid"
expSpace <- createExpSpace(attPerturb = attPerturb,
                           attPerturbSamp = attPerturbSamp,
                           attPerturbMin = attPerturbMin,
                           attPerturbMax = attPerturbMax,
                           attPerturbType = attPerturbType,
                           attHold = attHold,
                           attTargetsFile = NULL)

# perform stochastic simulation
simStochastic <- generateScenarios(reference = tank_obs,
                                   expSpace = expSpace,
                                   simLengthNyrs = 30,seedID = 1)

plotScenarios(simStochastic)

pause

#############################################################

# define the arguments of the systemModel, here tankWrapper
systemArgs <- list(roofArea = 205, 
                   nPeople = 1, 
                   tankVol = 2400, 
                   firstFlush = 2.0, 
                   write.file = FALSE)

metrics <- c("average daily deficit (L)", "reliability (fraction)", "volumetric reliability (fraction)")

systemModel = tankWrapper

clim = tank_obs

evaluate_system_metrics(sim=simScaling,
                        clim=clim,
                        systemModel=systemModel,
                        systemArgs=systemArgs,
                        metrics=metrics)
  
evaluate_system_metrics(sim=simStochastic,
                        clim=clim,
                        systemModel=systemModel,
                        systemArgs=systemArgs,
                        metrics=metrics)


#############################

# systemPerf <- runSystemModel(sim = sim,                     # simulation; the perturbed time series
#                              systemModel = systemModel,      # the system model function
#                              systemArgs = systemArgs,        # argument to the system model function
#                              metrics = metrics)              # selected performance metrics 
# 
# systemModel = tankWrapper



