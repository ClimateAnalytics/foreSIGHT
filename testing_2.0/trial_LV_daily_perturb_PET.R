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
modelSelection[["modelParameterBounds"]][["P"]][["alpha.SON"]] <- c(0.2, 0.9)
modelSelection[["modelParameterBounds"]][["P"]][["alpha.DJF"]] <- c(0.2, 0.9)
modelSelection[["modelParameterBounds"]][["P"]][["alpha.MAM"]] <- c(0.2, 0.9)
modelSelection[["modelParameterBounds"]][["P"]][["alpha.JJA"]] <- c(0.2, 0.9)
modelSelection[["modelParameterBounds"]][["P"]][["sigma.SON"]] <- c(0.1, 5)
modelSelection[["modelParameterBounds"]][["P"]][["sigma.DJF"]] <- c(0.1, 5)
modelSelection[["modelParameterBounds"]][["P"]][["sigma.MAM"]] <- c(0.1, 5)
modelSelection[["modelParameterBounds"]][["P"]][["sigma.JJA"]] <- c(0.1, 5)
modelSelection[["modelParameterBounds"]][["P"]][["mu.SON"]] <- c(-5, 1)
modelSelection[["modelParameterBounds"]][["P"]][["mu.DJF"]] <- c(-5, 1)
modelSelection[["modelParameterBounds"]][["P"]][["mu.MAM"]] <- c(-5, 1)
modelSelection[["modelParameterBounds"]][["P"]][["mu.JJA"]] <- c(-5, 1)
modelSelection[["modelParameterBounds"]][["P"]][["lambda.SON"]] <- c(1, 3)
modelSelection[["modelParameterBounds"]][["P"]][["lambda.DJF"]] <- c(1, 3)
modelSelection[["modelParameterBounds"]][["P"]][["lambda.MAM"]] <- c(1, 3)
modelSelection[["modelParameterBounds"]][["P"]][["lambda.JJA"]] <- c(1, 3)

modelSelection[["penaltyAttributes"]] <- c("P_day_all_tot_m",'P_day_all_P99','P_day_all_nWet','P_day_all_avgDSD_m')
#modelSelection[["penaltyAttributes"]] <- c("P_day_all_tot_m",'P_day_all_xP99overPave','P_day_all_nWet','P_day_all_avgDSD_m')
modelSelection[["penaltyWeights"]] <- c(2,2,2,2)

pars = c(0.6204190,0.5958725,0.6663199,0.4945255,1.3575024,2.7440208,
         1.4322516,1.5119136,-0.4585758,-3.1848170,-0.5564916,0.1969416,
         2.1777857,1.6794270,2.5668475,1.8902215)

modelSelection$optimisationArguments = list()
modelSelection$optimisationArguments$suggestions = pars
modelSelection$optimisationArguments$nMultiStart = 1

modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
write(modelSelectionJSON, file = controlFile)

####################

attPerturb = c("P_day_all_tot_m")
# attHold = c('P_day_all_P99','P_day_all_nWet','P_day_all_avgWSD_m',
#             "P_day_SON_tot_m",'P_day_SON_P99','P_day_SON_nWet','P_day_SON_avgWSD_m',
#             'P_day_DJF_tot_m','P_day_DJF_P99','P_day_DJF_nWet','P_day_DJF_avgWSD_m',
#             'P_day_MAM_tot_m','P_day_MAM_P99','P_day_MAM_nWet','P_day_MAM_avgWSD_m',
#             'P_day_JJA_tot_m','P_day_JJA_P99','P_day_JJA_nWet','P_day_JJA_avgWSD_m')

# attHold = c('P_day_DJF_fracTot','P_day_MAM_fracTot','P_day_JJA_fracTot','P_day_SON_fracTot',
#            'P_day_all_nWet_m','P_day_DJF_fracNwet','P_day_MAM_fracNwet','P_day_JJA_fracNwet','P_day_SON_fracNwet',
#            'P_day_all_xP99overPave','P_day_DJF_fracxP99overPave','P_day_MAM_fracxP99overPave','P_day_JJA_fracxP99overPave','P_day_SON_fracxP99overPave',
#            'P_day_all_avgDSD_m','P_day_DJF_avgDSD_m','P_day_MAM_avgDSD_m','P_day_JJA_avgDSD_m','P_day_SON_avgDSD_m')

# attPerturbType = "regGrid"
# attPerturbSamp = c(1)
# attPerturbMin = c(1)
# attPerturbMax = c(1)

attHold = c('P_day_all_P99','P_day_all_nWet','P_day_all_avgDSD_m')
#attHold = c('P_day_all_xP99overPave','P_day_all_nWet','P_day_all_avgDSD_m')

attPerturbType = "regGrid"
attPerturbSamp = c(5)
attPerturbMin = c(0.8)
attPerturbMax = c(1.2)
# attPerturbSamp = c(1)
# attPerturbMin = c(0.8)
# attPerturbMax = c(0.8)

# create the exposure space
expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType,
                          attHold = attHold)

attsAll = c(attPerturb,attHold)
for (att in attsAll){
  for (seas in c('DJF','MAM','JJA','SON')){
    att.seas = gsub('all',seas,att)
    expSpace$targetMat[att.seas] = expSpace$targetMat[att]
    expSpace$attPerturb = c(expSpace$attPerturb,att.seas)
  }
}
#attsAll = c(expSpace$attPerturb,expSpace$attHold)
#expSpace$targetMat = expSpace$targetMat[,attsAll]

####################

sim_seas = generateScenarios(reference = clim,
                                   expSpace = expSpace,
                                   controlFile = controlFile,
                                   seedID = 1)

plotScenarios(sim_seas)

attSel = c('P_day_all_tot_m','P_day_all_xP99overPave','P_day_all_nWet','P_day_all_avgWSD_m',
              "P_day_SON_tot_m",'P_day_SON_xP99overPave','P_day_SON_nWet','P_day_SON_avgDSD_m',
              'P_day_DJF_tot_m','P_day_DJF_xP99overPave','P_day_DJF_nWet','P_day_DJF_avgDSD_m',
              'P_day_MAM_tot_m','P_day_MAM_xP99overPave','P_day_MAM_nWet','P_day_MAM_avgDSD_m',
              'P_day_JJA_tot_m','P_day_JJA_xP99overPave','P_day_JJA_nWet','P_day_JJA_avgDSD_m')
clim_sim = clim; clim_sim$P = sim_seas[[1]][[1]]$P$sim
attsObs = calculateAttributes(clim,attSel=attSel)
attsSim = calculateAttributes(clim_sim,attSel=attSel)
attsSim/attsObs
pause


####################
