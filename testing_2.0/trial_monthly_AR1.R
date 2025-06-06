rm(list=ls())

devtools::load_all()

tank_obs = convert_climYMD_POSIXct(tank_obs)

######################
# generate monthly rainfall 

clim = tank_obs

P.zoo = zoo(clim$P,clim$times)
P.zoo.mon = aggregate(P.zoo,zoo::as.yearmon,sum)
times.mon = time(P.zoo.mon)
clim_mon = list(times=as.POSIXct(times.mon,tz = 'UTC'),
                P=coredata(P.zoo.mon))

####################
# setup seasonal monAR1 model

modelSelection = list()
modelSelection$modelType = list()
modelSelection$modelType$P = "monAR1"
modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "seas"
modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")
write(modelSelectionJSON, file = controlFile)

####################

attPerturb = c("P_month_DJF_tot_m")
attHold = c("P_month_DJF_P50","P_month_DJF_P90","P_month_DJF_cor",
            "P_month_MAM_tot_m","P_month_MAM_P50","P_month_MAM_P90","P_month_MAM_cor",
            "P_month_JJA_tot_m","P_month_JJA_P50","P_month_JJA_P90","P_month_JJA_cor",
            "P_month_SON_tot_m","P_month_SON_P50","P_month_SON_P90","P_month_SON_cor")

attPerturbType = "regGrid"
#attPerturbSamp = c(2)
#attPerturbMin = c(0.9)
#attPerturbMax = c(1)

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

sim_seas = generateScenarios(reference = clim_mon,
                             expSpace = expSpace,
                             controlFile = controlFile,
                             seedID=1)

plotScenarios(sim_seas)


####################
# setup harmonic monAR1 model

# modelSelection = list()
# modelSelection$modelType = list()
# modelSelection$modelType$P = "monAR1"
# modelSelection$modelParameterVariation$P = "har"
# modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
# controlFile = paste0(tempdir(), "\\eg_controlFile.json")
# write(modelSelectionJSON, file = controlFile)
# 
# sim_har = generateScenarios(reference = clim_mon,
#                              expSpace = expSpace,
#                              controlFile = controlFile,
#                              seedID=1)
# 
# plotScenarios(sim_har)

