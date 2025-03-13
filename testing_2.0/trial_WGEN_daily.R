rm(list=ls())

devtools::load_all()

clim = convert_climYMD_POSIXct(tank_obs)


####################

modelSelection = list()
modelSelection$modelType = list()
modelSelection$modelType$P = "wgenDM"
#modelSelection$modelType$P = "latent"
modelSelection$modelParameterVariation = list()
#modelSelection$modelParameterVariation$P = "annual"
modelSelection$modelParameterVariation$P = "ann"
#modelSelection$modelParameterVariation$P = "annDelta"
modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")

write(modelSelectionJSON, file = controlFile)

####################

attPerturb = c("P_day_all_tot_m")
attHold = c('P_day_all_P99','P_day_all_nWet','P_day_all_avgWSD')

# attPerturb = c("P_day_JJA_tot_m")
# attHold = c('P_day_JJA_P99','P_day_JJA_nWet','P_day_JJA_avgWSD')

attPerturbType = "regGrid"
attPerturbSamp = c(1)
attPerturbMin = c(1.2)
attPerturbMax = c(1.2)

# create the exposure space
expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType,
                          attHold = attHold)


####################

# sim = generateScenarios(reference = clim,
#                         expSpace = expSpace,
#                         controlFile = controlFile,
#                         seedID=1)
# 
# plotScenarios(sim)

####################

#SEAS MODEL

modelSelection$modelParameterVariation$P = "seas"
modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)

# modelSelection[["optimisationArguments"]] <- list()
# modelSelection[["optimisationArguments"]][["nMultiStart"]] <- 10

modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
write(modelSelectionJSON, file = controlFile)

####################

attPerturb = c("P_day_all_tot_m")
attHold = c('P_day_all_P99','P_day_all_nWet','P_day_all_avgWSD',
            "P_day_SON_tot_m",'P_day_SON_P99','P_day_SON_nWet','P_day_SON_avgWSD',
            'P_day_DJF_tot_m','P_day_DJF_P99','P_day_DJF_nWet','P_day_DJF_avgWSD',
            'P_day_MAM_tot_m','P_day_MAM_P99','P_day_MAM_nWet','P_day_MAM_avgWSD',
            'P_day_JJA_tot_m','P_day_JJA_P99','P_day_JJA_nWet','P_day_JJA_avgWSD')

attPerturbType = "regGrid"
attPerturbSamp = c(1)
attPerturbMin = c(1)
attPerturbMax = c(1)

# create the exposure space
expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType,
                          attHold = attHold)

####################

# sim_seas = generateScenarios(reference = clim,
#                                    expSpace = expSpace,
#                                    controlFile = controlFile,
#                                    seedID = 1)
# 
# plotScenarios(sim_seas)

####################

#HAR MODEL

modelSelection$modelParameterVariation$P = "har"
modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
write(modelSelectionJSON, file = controlFile)

sim_har = generateScenarios(reference = clim,
                             expSpace = expSpace,
                             controlFile = controlFile,
                             seedID = 1)

plotScenarios(sim_har)