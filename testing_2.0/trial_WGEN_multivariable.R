rm(list=ls())

devtools::load_all()

clim = convert_climYMD_POSIXct(tank_obs)
  
######################################################################

# # Selected attributes
# attPerturb <- c("P_day_all_tot_m", "Temp_day_all_avg_m")
# attHold <- c("P_day_all_R10_m", "P_day_DJF_tot_m","Temp_day_all_rng_m", "Temp_day_DJF_avg_m")
# 
# # Sampling bounds and strategy
# attPerturbType = "regGrid"
# attPerturbSamp = c(2, 2)
# attPerturbMin = c(0.8,-0.5)
# attPerturbMax = c(1.2,0.5)

attPerturb <- c("Temp_day_all_avg_m")
#attHold <- c("Temp_day_all_rng_m", "Temp_day_DJF_avg_m")
attHold <- c("Temp_day_all_rng_m", "Temp_day_all_cor")

attAll = c(attPerturb,attHold)
calculateAttributes(clim,attAll)

# Sampling bounds and strategy
attPerturbType = "regGrid"
attPerturbSamp = c(2)
attPerturbMin = c(0.8)
attPerturbMax = c(1.2)

# Creating the exposure space
expSpace <- createExpSpace(attPerturb = attPerturb, 
                           attPerturbSamp = attPerturbSamp, 
                           attPerturbMin = attPerturbMin, 
                           attPerturbMax = attPerturbMax,
                           attPerturbType = attPerturbType, 
                           attHold = attHold)

modelSelection = list()
modelSelection$modelType = list()
modelSelection$modelType$Temp = "wgenDMtemp"
modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$Temp = "ann"
modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")

write(modelSelectionJSON, file = controlFile)

sim <- generateScenarios(reference = clim,   # reference time series 
                         expSpace = expSpace,    # exposure space
                         numReplicates = 3,      # number of replicates
                         controlFile = controlFile)


plotScenarios(sim)

####################

#modelSelection$modelParameterVariation$Temp = "har"
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

