rm(list=ls())

devtools::load_all()
clim = convert_climYMD_POSIXct(tank_obs)
#clim = tank_obs

mvFunc_cor = function(data.1,data.2){
  return(cor(data.1,data.2,use='pairwise.complete.obs'))
}

mvFunc_meanWetDay = function(data.1,data.2){
  return(mean(data.1[data.2>0]))
}

mvFunc_sdWetDay = function(data.1,data.2){
  return(sd(data.1[data.2>0]))
}

mvFunc_meanDryDay = function(data.1,data.2){
  return(mean(data.1[data.2==0]))
}

mvFunc_sdDryDay = function(data.1,data.2){
  return(sd(data.1[data.2==0]))
}

attPerturb = c("P_day_all_tot_m")
attHold = c('P_day_all_P99','P_day_all_nWet','P_day_all_avgWSD_m',
            "P_day_SON_tot_m",'P_day_SON_P99','P_day_SON_nWet','P_day_SON_avgWSD_m',
            'P_day_DJF_tot_m','P_day_DJF_P99','P_day_DJF_nWet','P_day_DJF_avgWSD_m',
            'P_day_MAM_tot_m','P_day_MAM_P99','P_day_MAM_nWet','P_day_MAM_avgWSD_m',
            'P_day_JJA_tot_m','P_day_JJA_P99','P_day_JJA_nWet','P_day_JJA_avgWSD_m',
            "Temp_day_all_cor","mv.Temp.P_day_all_meanWetDay","mv.Temp.P_day_all_sdWetDay",
            "mv.Temp.P_day_all_meanDryDay","mv.Temp.P_day_all_sdDryDay")
attSel = c(attPerturb,attHold)
calculateAttributes(clim,attSel=attSel)

# attPerturb = c("Temp_day_all_avg_m","P_day_all_tot_m")
# attHold = c('P_day_all_P99','P_day_all_nWet','P_day_all_avgWSD_m',
#             "Temp_day_all_rng_m", "Temp_day_all_cor",
#             'mv.Temp.P_day_all_cor','mv.Temp.P_day_all_meanWD')

# attPerturb = c("Temp_day_all_avg_m","P_day_all_tot_m")
# attHold = c('P_day_all_P99','P_day_all_nWet','P_day_all_avgWSD_m',
#             "Temp_day_all_rng_m", "Temp_day_all_cor",
#             'mv.Temp.P_day_all_cor','mv.Temp.P_day_all_meanWD')

# attPerturb = c("P_day_all_tot_m")
# attHold = c('P_day_all_P99','P_day_all_nWet','P_day_all_avgWSD_m',
#             "Temp_day_all_cor",
#             "mv.Temp.P_day_all_meanWetDay","mv.Temp.P_day_all_sdWetDay",
#             "mv.Temp.P_day_all_meanDryDay","mv.Temp.P_day_all_sdDryDay")

# attPerturb = c("P_day_all_tot_m")
# attHold = c('P_day_all_P99','P_day_all_nWet','P_day_all_avgWSD_m')

# devtools::load_all('C:/Users/a1065639/Work/foreSIGHT1.2/foreSIGHT/')
# clim = tank_obs 
# attPerturb = c("P_ann_tot_m","Temp_ann_avg_m")
# attHold = c('P_ann_P99','P_ann_nWet','P_ann_avgWSD_m',
#             "Temp_ann_rng_m", "Temp_ann_cor")

######################################################################

# 
# 
# 
# 
# 
# attPerturbType = "regGrid"
# # attPerturbSamp = c(1,1)
# # attPerturbMin = c(1,1)
# # attPerturbMax = c(1,1)
# attPerturbSamp = c(1)
# attPerturbMin = c(1)
# attPerturbMax = c(1)
# 
# # Creating the exposure space
# expSpace <- createExpSpace(attPerturb = attPerturb, 
#                            attPerturbSamp = attPerturbSamp, 
#                            attPerturbMin = attPerturbMin, 
#                            attPerturbMax = attPerturbMax,
#                            attPerturbType = attPerturbType, 
#                            attHold = attHold)
# 
# ######################################################################
# 
# modelSelection = list()
# modelSelection$modelType = list()
# modelSelection$modelParameterVariation = list()
# 
# modelSelection$modelType$P = "wgen"
# modelSelection$modelType$Temp = "wgenLM"
# modelSelection$modelParameterVariation$P = "ann"
# modelSelection$modelParameterVariation$Temp = "annWD"
# 
# modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
# controlFile = paste0(tempdir(), "\\eg_controlFile.json")
# 
# write(modelSelectionJSON, file = controlFile)
# 
# sim <- generateScenarios(reference = clim,   # reference time series 
#                          expSpace = expSpace,    # exposure space
#                          numReplicates = 1,      # number of replicates
#                          controlFile = controlFile,simLengthNyrs = 30)
# 
# 
# plotScenarios(sim)

####################

attPerturb = c("P_day_all_tot_m")
attHold = c('P_day_all_P99','P_day_all_nWet','P_day_all_avgWSD_m',
            "P_day_SON_tot_m",'P_day_SON_P99','P_day_SON_nWet','P_day_SON_avgWSD_m',
            'P_day_DJF_tot_m','P_day_DJF_P99','P_day_DJF_nWet','P_day_DJF_avgWSD_m',
            'P_day_MAM_tot_m','P_day_MAM_P99','P_day_MAM_nWet','P_day_MAM_avgWSD_m',
            'P_day_JJA_tot_m','P_day_JJA_P99','P_day_JJA_nWet','P_day_JJA_avgWSD_m',
            "Temp_day_all_cor","mv.Temp.P_day_all_meanWetDay","mv.Temp.P_day_all_sdWetDay",
            "mv.Temp.P_day_all_meanDryDay","mv.Temp.P_day_all_sdDryDay")
            # "Temp_day_all_cor","mv.Temp.P_day_all_meanWetDay","mv.Temp.P_day_all_sdWetDay",
            # "mv.Temp.P_day_all_meanDryDay","mv.Temp.P_day_all_sdDryDay",
            # "Temp_day_SON_cor","mv.Temp.P_day_SON_meanWetDay","mv.Temp.P_day_SON_sdWetDay",
            # "mv.Temp.P_day_SON_meanDryDay","mv.Temp.P_day_SON_sdDryDay",
            # "Temp_day_DJF_cor","mv.Temp.P_day_DJF_meanWetDay","mv.Temp.P_day_DJF_sdWetDay",
            # "mv.Temp.P_day_DJF_meanDryDay","mv.Temp.P_day_DJF_sdDryDay",
            # "Temp_day_MAM_cor","mv.Temp.P_day_MAM_meanWetDay","mv.Temp.P_day_MAM_sdWetDay",
            # "mv.Temp.P_day_MAM_meanDryDay","mv.Temp.P_day_MAM_sdDryDay",
            # "Temp_day_JJA_cor","mv.Temp.P_day_JJA_meanWetDay","mv.Temp.P_day_JJA_sdWetDay",
            # "mv.Temp.P_day_JJA_meanDryDay","mv.Temp.P_day_JJA_sdDryDay")


modelSelection = list()
modelSelection$modelType = list()
modelSelection$modelParameterVariation = list()
#modelSelection$modelType$P = "wgen"
modelSelection$modelType$P = "latent"
modelSelection$modelParameterVariation$P = "seas"
modelSelection$modelType$Temp = "wgenLM"
modelSelection$modelParameterVariation$Temp = "annWD"
modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")

write(modelSelectionJSON, file = controlFile)

# Sampling bounds and strategy
attPerturbType = "regGrid"
attPerturbSamp = c(1)
attPerturbMin = c(1)
attPerturbMax = c(1)

# Creating the exposure space
expSpace <- createExpSpace(attPerturb = attPerturb, 
                           attPerturbSamp = attPerturbSamp, 
                           attPerturbMin = attPerturbMin, 
                           attPerturbMax = attPerturbMax,
                           attPerturbType = attPerturbType, 
                           attHold = attHold)

sim_seas <- generateScenarios(reference = clim,   # reference time series 
                            expSpace = expSpace,    # exposure space
                            numReplicates = 1,      # number of replicates
                            controlFile = controlFile,
                            simLengthNyrs = 100
                            )

plotScenarios(sim_seas)

pause

####################


modelSelection$modelParameterVariation$Temp = "seas"
modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")
write(modelSelectionJSON, file = controlFile)

attPerturb <- c("Temp_day_SON_avg_m")
attHold <- c("Temp_day_SON_rng_m", "Temp_day_SON_cor",
             "Temp_day_DJF_avg_m","Temp_day_DJF_rng_m", "Temp_day_DJF_cor",
             "Temp_day_MAM_avg_m","Temp_day_MAM_rng_m", "Temp_day_MAM_cor",
             "Temp_day_JJA_avg_m","Temp_day_JJA_rng_m", "Temp_day_JJA_cor")

attAll = c(attPerturb,attHold)
calculateAttributes(clim,attAll)

# Sampling bounds and strategy
attPerturbType = "regGrid"
attPerturbSamp = c(1)
attPerturbMin = c(0)
attPerturbMax = c(0)

# Creating the exposure space
expSpace <- createExpSpace(attPerturb = attPerturb, 
                           attPerturbSamp = attPerturbSamp, 
                           attPerturbMin = attPerturbMin, 
                           attPerturbMax = attPerturbMax,
                           attPerturbType = attPerturbType, 
                           attHold = attHold)

sim_har<- generateScenarios(reference = clim,   # reference time series 
                         expSpace = expSpace,    # exposure space
                         numReplicates = 3,      # number of replicates
                         controlFile = controlFile)


plotScenarios(sim_har)

####################

