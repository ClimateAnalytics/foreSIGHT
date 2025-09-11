rm(list=ls())

devtools::load_all()

#tank_obs = convert_climYMD_POSIXct(tank_obs)

######################
# generate monthly rainfall 

clim = tank_obs

P.zoo = zoo::zoo(clim$P,clim$times)
P.zoo.mon = aggregate(P.zoo,zoo::as.yearmon,sum)
times.mon = time(P.zoo.mon)
clim_mon = list(times=as.POSIXct(times.mon,tz = 'UTC'),
                P=zoo::coredata(P.zoo.mon))

####################
# setup seasonal monAR1 model

modelSelection = list()
modelSelection$modelType = list()
modelSelection$modelType$P = "monAR1"
modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "seas1"

# modelSelection[["penaltyAttributes"]] <- c("P_month_all_tot_m","P_month_all_P50","P_month_all_P90","P_month_all_cor")
# modelSelection[["penaltyWeights"]] = rep(3,4)

modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")
write(modelSelectionJSON, file = controlFile)

####################

func_skew = function(data){
  return(moments::skewness(data,na.rm = T))
}

####################

# attPerturb = c("P_month_DJF_tot_m")
# # attHold = c("P_month_DJF_P50","P_month_DJF_P90","P_month_DJF_cor",
# #             "P_month_MAM_tot_m","P_month_MAM_P50","P_month_MAM_P90","P_month_MAM_cor",
# #             "P_month_JJA_tot_m","P_month_JJA_P50","P_month_JJA_P90","P_month_JJA_cor",
# #             "P_month_SON_tot_m","P_month_SON_P50","P_month_SON_P90","P_month_SON_cor")
# 
# attHold = c("P_month_DJF_P50","P_month_DJF_P90",
#             "P_month_MAM_tot_m","P_month_MAM_P50","P_month_MAM_P90",
#             "P_month_JJA_tot_m","P_month_JJA_P50","P_month_JJA_P90",
#             "P_month_SON_tot_m","P_month_SON_P50","P_month_SON_P90",
#             "P_month_all_cor")


# attPerturb = c("P_month_all_tot_m")
# attHold = c("P_month_all_cv","P_month_all_skew","P_month_all_cor",
#             "P_month_DJF_tot_m","P_month_DJF_cv","P_month_DJF_skew","P_month_DJF_cor",
#             "P_month_MAM_tot_m","P_month_MAM_cv","P_month_MAM_skew","P_month_MAM_cor",
#             "P_month_JJA_tot_m","P_month_JJA_cv","P_month_JJA_skew","P_month_JJA_cor",
#             "P_month_SON_tot_m","P_month_SON_cv","P_month_SON_skew","P_month_SON_cor")


attPerturb = c("P_month_DJF_tot_m")
attHold = c("P_month_DJF_cv","P_month_DJF_skew",
            "P_month_MAM_tot_m","P_month_MAM_cv","P_month_MAM_skew",
            "P_month_JJA_tot_m","P_month_JJA_cv","P_month_JJA_skew",
            "P_month_SON_tot_m","P_month_SON_cv","P_month_SON_skew",
            "P_month_all_cor")

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
                             seedID=1,simLengthNyrs = 100)

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

