rm(list=ls())

devtools::load_all()

####################

times = as.POSIXct(paste0(tank_obs$year,'/',
                          tank_obs$month,'/',
                          tank_obs$day),
                   tz = 'UTC')

clim = list(times=times,
            P=tank_obs$P)

####################

calculateAttributes(clim,'P_day_SON_tot')

####################

modelSelection = list()
modelSelection$modelType = list()
modelSelection$modelType$P = "wgenDM"
modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "seas"
modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")

write(modelSelectionJSON, file = controlFile)

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

sim_seas = generateScenarios(reference = clim,
                                   expSpace = expSpace,
                                   controlFile = controlFile,
                                   seedID = 1,simLengthNyrs = 100)

plotScenarios(sim_seas)

