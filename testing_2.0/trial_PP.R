rm(list=ls())

devtools::load_all()

clim = convert_climYMD_POSIXct(tank_obs)


####################

modelSelection = list()
modelSelection$modelType = list()
#modelSelection$modelType$P = "wgenDM"
modelSelection$modelType$P = "latent"
modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "seas"
modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")

write(modelSelectionJSON, file = controlFile)

####################

attPerturb = c("P_day_all_tot_m")
attHold = c('P_day_all_xP99overPave','P_day_all_nWet','P_day_all_avgWSD',
            "P_day_SON_tot_m",'P_day_SON_xP99overPave','P_day_SON_nWet','P_day_SON_avgWSD',
            'P_day_DJF_tot_m','P_day_DJF_xP99overPave','P_day_DJF_nWet','P_day_DJF_avgWSD',
            'P_day_MAM_tot_m','P_day_MAM_xP99overPave','P_day_MAM_nWet','P_day_MAM_avgWSD',
            'P_day_JJA_tot_m','P_day_JJA_xP99overPave','P_day_JJA_nWet','P_day_JJA_avgWSD')

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

clim_sim = list(times=sim_seas$simDates,P=sim_seas$Rep1$Target1$P$sim)

calculateAttributes(clim,'P_day_all_P99.9')
calculateAttributes(clim_sim,'P_day_all_P99.9')

calculateAttributes(clim,'P_day_all_P99')
calculateAttributes(clim_sim,'P_day_all_P99')

pause

####################

modelSelection$modelType$P = "latent"
modelSelection$postProcessing = list(P=list())
modelSelection$postProcessing$P$types = c('scaleExtremesAll')

modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
write(modelSelectionJSON, file = controlFile)

attPerturb = c("P_day_all_tot_m")
attHold = c('P_day_all_xP99overPave','P_day_all_nWet','P_day_all_avgWSD',
            "P_day_SON_tot_m",'P_day_SON_xP99overPave','P_day_SON_nWet','P_day_SON_avgWSD',
            'P_day_DJF_tot_m','P_day_DJF_xP99overPave','P_day_DJF_nWet','P_day_DJF_avgWSD',
            'P_day_MAM_tot_m','P_day_MAM_xP99overPave','P_day_MAM_nWet','P_day_MAM_avgWSD',
            'P_day_JJA_tot_m','P_day_JJA_xP99overPave','P_day_JJA_nWet','P_day_JJA_avgWSD')

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


sim_seas_scaleExtremesAll = generateScenarios(reference = clim,
                             expSpace = expSpace,
                             controlFile = controlFile,
                             seedID = 1)

plotScenarios(sim_seas_scaleExtremesAll)

clim_sim_scaleAll = list(times=sim_seas_scaleExtremesAll$simDates,P=sim_seas_scaleExtremesAll$Rep1$Target1$P$sim)

calculateAttributes(clim,'P_day_all_P99.9')
calculateAttributes(clim_sim_scaleAll,'P_day_all_P99.9')

calculateAttributes(clim,'P_day_all_P99')
calculateAttributes(clim_sim_scaleAll,'P_day_all_P99')

calculateAttributes(clim,'P_day_DJF_P99.9')
calculateAttributes(clim_sim_scaleAll,'P_day_DJF_P99.9')


browser()

####################

modelSelection$modelType$P = "latent"
modelSelection$postProcessing = list(P=list())
modelSelection$postProcessing$P$types = c('scaleExtremesSeas')

modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
write(modelSelectionJSON, file = controlFile)

attPerturb = c("P_day_all_tot_m")
attHold = c('P_day_all_xP99overPave','P_day_all_nWet','P_day_all_avgWSD',
            "P_day_SON_tot_m",'P_day_SON_xP99overPave','P_day_SON_nWet','P_day_SON_avgWSD',
            'P_day_DJF_tot_m','P_day_DJF_xP99overPave','P_day_DJF_nWet','P_day_DJF_avgWSD',
            'P_day_MAM_tot_m','P_day_MAM_xP99overPave','P_day_MAM_nWet','P_day_MAM_avgWSD',
            'P_day_JJA_tot_m','P_day_JJA_xP99overPave','P_day_JJA_nWet','P_day_JJA_avgWSD')

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


sim_seas_scaleExtremesSeas = generateScenarios(reference = clim,
                                              expSpace = expSpace,
                                              controlFile = controlFile,
                                              seedID = 1)

plotScenarios(sim_seas_scaleExtremesSeas)

clim_sim_scaleSeas = list(times=sim_seas_scaleExtremesSeas$simDates,P=sim_seas_scaleExtremesSeas$Rep1$Target1$P$sim)

calculateAttributes(clim,'P_day_all_P99.9')
calculateAttributes(clim_sim_scaleSeas,'P_day_all_P99.9')

calculateAttributes(clim,'P_day_all_P99')
calculateAttributes(clim_sim_scaleSeas,'P_day_all_P99')

calculateAttributes(clim,'P_day_DJF_P99.9')
calculateAttributes(clim_sim_scaleSeas,'P_day_DJF_P99.9')

####################

modelSelection$postProcessing = list()
modelSelection$postProcessing$P = list()
modelSelection$postProcessing$P$types = c('annVar')

modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
write(modelSelectionJSON, file = controlFile)

attPerturb = c("P_day_all_tot_m")
attHold = c('P_day_all_xP99overPave','P_day_all_nWet','P_day_all_avgWSD',
            "P_day_SON_tot_m",'P_day_SON_xP99overPave','P_day_SON_nWet','P_day_SON_avgWSD',
            'P_day_DJF_tot_m','P_day_DJF_xP99overPave','P_day_DJF_nWet','P_day_DJF_avgWSD',
            'P_day_MAM_tot_m','P_day_MAM_xP99overPave','P_day_MAM_nWet','P_day_MAM_avgWSD',
            'P_day_JJA_tot_m','P_day_JJA_xP99overPave','P_day_JJA_nWet','P_day_JJA_avgWSD',
            'P_year_all_sd')
# attHold = c('P_day_all_P99','P_day_all_nWet','P_day_all_avgWSD',
#             "P_day_SON_tot_m",'P_day_SON_P99','P_day_SON_nWet','P_day_SON_avgWSD',
#             'P_day_DJF_tot_m','P_day_DJF_P99','P_day_DJF_nWet','P_day_DJF_avgWSD',
#             'P_day_MAM_tot_m','P_day_MAM_P99','P_day_MAM_nWet','P_day_MAM_avgWSD',
#             'P_day_JJA_tot_m','P_day_JJA_P99','P_day_JJA_nWet','P_day_JJA_avgWSD',
#             'P_year_all_sd')

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

func_sd = function(data) sd(data)

sim_seas_annSD = generateScenarios(reference = clim,
                                   expSpace = expSpace,
                                   controlFile = controlFile,
                                   seedID = 1)

plotScenarios(sim_seas_annSD)

####################

