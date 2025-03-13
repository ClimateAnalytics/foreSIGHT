rm(list=ls())

devtools::load_all()

clim = convert_climYMD_POSIXct(tank_obs)


####################

modelSelection = list()
modelSelection$modelType = list()
modelSelection$modelType$P = "LV"
#modelSelection$modelType$P = "latent"
modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "ann"
modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")

write(modelSelectionJSON, file = controlFile)

####################

# #attPerturb = c("P_day_all_tot_m")
# #attHold = c('P_day_all_P99','P_day_all_nWet','P_day_all_cor')
# 
# attPerturb = c("P_day_JJA_tot_m")
# attHold = c('P_day_JJA_P99','P_day_JJA_nWet','P_day_JJA_avgWSD_m')
# 
# # atts = calculateAttributes(clim,attSel=c(attPerturb,attHold))
# #
# # pause
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


####################

# sim = generateScenarios(reference = clim,
#                         expSpace = expSpace,
#                         controlFile = controlFile,
#                         seedID=1)
# 
# plotScenarios(sim)
# 
# pause

####################

#SEAS MODEL

modelSelection$modelParameterVariation$P = "seas"

# modelSelection[["optimisationArguments"]] <- list()
# modelSelection[["optimisationArguments"]][["nMultiStart"]] <- 10

modelSelection[["modelParameterBounds"]] <- list()
modelSelection[["modelParameterBounds"]][["P"]] <- list()
modelSelection[["modelParameterBounds"]][["P"]][["lambda.SON"]] <- c(0.5, 4)
modelSelection[["modelParameterBounds"]][["P"]][["lambda.DJF"]] <- c(0.5, 4)
modelSelection[["modelParameterBounds"]][["P"]][["lambda.MAM"]] <- c(0.5, 4)
modelSelection[["modelParameterBounds"]][["P"]][["lambda.JJA"]] <- c(0.5, 4)
modelSelection[["modelParameterBounds"]][["P"]][["mu.SON"]] <- c(-15, 1)
modelSelection[["modelParameterBounds"]][["P"]][["mu.DJF"]] <- c(-15, 1)
modelSelection[["modelParameterBounds"]][["P"]][["mu.MAM"]] <- c(-15, 1)
modelSelection[["modelParameterBounds"]][["P"]][["mu.JJA"]] <- c(-15, 1)

modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
write(modelSelectionJSON, file = controlFile)

####################

attPerturb = c("P_day_all_tot_m")
attHold = c('P_day_all_P99','P_day_all_nWet','P_day_all_avgWSD_m',
            "P_day_SON_tot_m",'P_day_SON_P99','P_day_SON_nWet','P_day_SON_avgWSD_m',
            'P_day_DJF_tot_m','P_day_DJF_P99','P_day_DJF_nWet','P_day_DJF_avgWSD_m',
            'P_day_MAM_tot_m','P_day_MAM_P99','P_day_MAM_nWet','P_day_MAM_avgWSD_m',
            'P_day_JJA_tot_m','P_day_JJA_P99','P_day_JJA_nWet','P_day_JJA_avgWSD_m')

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

sim_seas = generateScenarios(reference = clim,
                                   expSpace = expSpace,
                                   controlFile = controlFile,
                                   seedID = 1)

plotScenarios(sim_seas)

save.image(file='testing_2.0/sim_seas_LV.RData')

####################
