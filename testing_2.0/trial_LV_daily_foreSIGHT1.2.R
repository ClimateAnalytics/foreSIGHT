rm(list=ls())

devtools::load_all()

# clim = convert_climYMD_POSIXct(tank_obs)


####################

modelSelection = list()
modelSelection$modelType = list()
#modelSelection$modelType$P = "LV"
modelSelection$modelType$P = "latent"
modelSelection$modelParameterVariation = list()
#modelSelection$modelParameterVariation$P = "ann"
modelSelection$modelParameterVariation$P = "annual"
modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")

write(modelSelectionJSON, file = controlFile)

####################

# attPerturb = c("P_JJA_tot_m")
# attHold = c('P_JJA_P99','P_JJA_nWet','P_JJA_avgWSD_m')
# 
# attPerturbType = "regGrid"
# attPerturbSamp = c(1)
# attPerturbMin = c(1)
# attPerturbMax = c(1)
# 
# # create the exposure space
# expSpace = createExpSpace(attPerturb = attPerturb,
#                           attPerturbSamp = attPerturbSamp,
#                           attPerturbMin = attPerturbMin,
#                           attPerturbMax = attPerturbMax,
#                           attPerturbType = attPerturbType,
#                           attHold = attHold)
# 
# 
# ####################
# 
# sim = generateScenarios(reference = tank_obs,
#                         expSpace = expSpace,
#                         controlFile = controlFile,
#                         seedID=1)
# 
# plotScenarios(sim)
# 
# pause

####################

#SEAS MODEL

#modelSelection$modelParameterVariation$P = "seas"
modelSelection$modelParameterVariation$P = "seasonal"

modelSelection[["modelParameterBounds"]] <- list()
modelSelection[["modelParameterBounds"]][["P"]] <- list()
modelSelection[["modelParameterBounds"]][["P"]][["lambda_1"]] <- c(0.5, 4)
modelSelection[["modelParameterBounds"]][["P"]][["lambda_2"]] <- c(0.5, 4)
modelSelection[["modelParameterBounds"]][["P"]][["lambda_3"]] <- c(0.5, 4)
modelSelection[["modelParameterBounds"]][["P"]][["lambda_4"]] <- c(0.5, 4)

# modelSelection[["optimisationArguments"]] <- list()
# modelSelection[["optimisationArguments"]][["nMultiStart"]] <- 10

modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
write(modelSelectionJSON, file = controlFile)

####################

attPerturb = c("P_ann_tot_m")
attHold = c('P_ann_P99','P_ann_nWet','P_ann_avgWSD_m',
            "P_SON_tot_m",'P_SON_P99','P_SON_nWet','P_SON_avgWSD_m',
            'P_DJF_tot_m','P_DJF_P99','P_DJF_nWet','P_DJF_avgWSD_m',
            'P_MAM_tot_m','P_MAM_P99','P_MAM_nWet','P_MAM_avgWSD_m',
            'P_JJA_tot_m','P_JJA_P99','P_JJA_nWet','P_JJA_avgWSD_m')

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

sim_seas = generateScenarios(reference = tank_obs,
                                   expSpace = expSpace,
                                   controlFile = controlFile,
                                   seedID = 1)

plotScenarios(sim_seas)

