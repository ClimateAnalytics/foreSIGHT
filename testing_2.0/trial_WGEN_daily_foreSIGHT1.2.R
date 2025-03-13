rm(list=ls())

devtools::load_all('C:/Users/a1065639/Work/foreSIGHT1.2/foreSIGHT/')



####################

modelSelection = list()
modelSelection$modelType = list()
modelSelection$modelType$P = "wgen"
modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "annual"
modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")

write(modelSelectionJSON, file = controlFile)

####################

attPerturb = c("P_ann_tot_m")
attHold = c('P_ann_P99','P_ann_nWet','P_ann_avgWSD')

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

# sim = generateScenarios(reference = tank_obs,
#                         expSpace = expSpace,
#                         controlFile = controlFile,
#                         seedID=1)
# 
# plotScenarios(sim)

####################

#SEAS MODEL

modelSelection$modelParameterVariation$P = "seasonal"
modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)

# modelSelection[["optimisationArguments"]] <- list()
# modelSelection[["optimisationArguments"]][["nMultiStart"]] <- 10

modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
write(modelSelectionJSON, file = controlFile)

####################

attPerturb = c("P_ann_tot_m")
attHold = c('P_ann_P99','P_ann_nWet','P_ann_avgWSD',
            "P_SON_tot_m",'P_SON_P99','P_SON_nWet','P_SON_avgWSD',
            'P_DJF_tot_m','P_DJF_P99','P_DJF_nWet','P_DJF_avgWSD',
            'P_MAM_tot_m','P_MAM_P99','P_MAM_nWet','P_MAM_avgWSD',
            'P_JJA_tot_m','P_JJA_P99','P_JJA_nWet','P_JJA_avgWSD')

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

# sim_seas = generateScenarios(reference = tank_obs,
#                                    expSpace = expSpace,
#                                    controlFile = controlFile,
#                                    seedID = 1)
# 
# plotScenarios(sim_seas)

####################

#HAR MODEL

modelSelection$modelParameterVariation$P = "harmonic"
modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
write(modelSelectionJSON, file = controlFile)

sim_har = generateScenarios(reference = tank_obs,
                            expSpace = expSpace,
                            controlFile = controlFile,
                            seedID = 1)

plotScenarios(sim_har)
