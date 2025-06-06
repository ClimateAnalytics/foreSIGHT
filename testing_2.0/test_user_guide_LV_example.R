devtools::load_all()

# create the exposure space
attPerturb <- c("P_day_all_tot_m", "P_day_all_P99")
attHold <- c("P_day_all_maxWSD_m", "P_day_all_nWet_m")
attPerturbType = "regGrid"
attPerturbSamp = c(3, 3)
attPerturbMin = c(0.8, 0.8)
attPerturbMax = c(1.2, 1.2)
expSpace <- createExpSpace(attPerturb = attPerturb, 
                           attPerturbSamp = attPerturbSamp, 
                           attPerturbMin = attPerturbMin, 
                           attPerturbMax = attPerturbMax,
                           attPerturbType = attPerturbType, 
                           attHold = attHold)

# specify the penalty settings in a list
controlFileList <- list()
controlFileList[["penaltyAttributes"]] <- c("P_day_all_tot_m")
controlFileList[["penaltyWeights"]] <- c(0.5)

# specify the alternate model selections
controlFileList[["modelType"]] <- list()
controlFileList[["modelType"]][["P"]] <- "latent"
controlFileList[["modelParameterVariation"]] <- list()
controlFileList[["modelParameterVariation"]][["P"]] <- "har"

# add user-specified bounds for model parameters
controlFileList[["modelParameterBounds"]] <- list()
controlFileList[["modelParameterBounds"]][["P"]] <- list()
controlFileList[["modelParameterBounds"]][["P"]][["mu.m"]] <- c(-5, 0)
controlFileList[["modelParameterBounds"]][["P"]][["alpha.m"]] <- c(0.35, 0.95)

# write the list into a JSON file
controlFileJSON <- jsonlite::toJSON(controlFileList, pretty = TRUE, auto_unbox = TRUE)
write(controlFileJSON, file = paste0(tempdir(), "controlFile.json"))

# generate scenarios
clim = convert_climYMD_POSIXct(tank_obs)

sim <- generateScenarios(reference = clim, expSpace = expSpace, 
                         controlFile = paste0(tempdir(), "controlFile.json"),
                         numReplicates = 10)



############

# attSel = c(attPerturb,attHold)

# attOther = c('P_day_all_tot_m','P_day_SON_tot_m','P_day_DJF_tot_m','P_day_MAM_tot_m','P_day_JJA_tot_m',
#              'P_day_all_P99','P_day_SON_P99','P_day_DJF_P99','P_day_MAM_P99','P_day_JJA_P99')

# attSel = unique(c(attPerturb,attHold,attOther))
# attSel = unique(c(attPerturb,attOther))

attSel = c('P_day_all_tot_m','P_day_SON_tot_m','P_day_DJF_tot_m','P_day_MAM_tot_m','P_day_JJA_tot_m',
           'P_day_all_P99','P_day_SON_P99','P_day_DJF_P99','P_day_MAM_P99','P_day_JJA_P99',
           'P_day_all_P99.9','P_day_SON_P99.9','P_day_DJF_P99.9','P_day_MAM_P99.9','P_day_JJA_P99.9',
           'P_day_all_nWet_m','P_day_SON_nWet_m','P_day_DJF_nWet_m','P_day_MAM_nWet_m','P_day_JJA_nWet_m',
           'P_day_all_maxWSD_m','P_day_SON_maxWSD_m','P_day_DJF_maxWSD_m','P_day_MAM_maxWSD_m','P_day_JJA_maxWSD_m')

P = calcPerformanceAttributes(clim=clim,sim=sim,attSel=attSel)

attSel='P_day_all_tot_m'
par(mfrow=c(5,5),mar=c(4,5,2,1))
for (att in names(P)){
  # plotPerformanceSpace(P, sim, metric=metric)
  plotPerformanceOAT(P, sim, metric=att,col='black',use_ggplot = F,attSel=attSel)
}

attSel='P_day_all_P99'
par(mfrow=c(5,5),mar=c(4,4,2,1))
for (att in names(P)){
  plotPerformanceOAT(P, sim, metric=att,col='black',use_ggplot = F,attSel=attSel)
}

pause

names(P) = paste0('.',names(P))

attSel='P_day_all_tot_m'
#par(mfrow=c(4,4),mar=c(4,4,2,1))
for (metric in names(P)){
  tar=sim$expSpace$targetMat
  att = strsplit(metric,'[.]')[[1]][2]
  if (att%in%names(tar)){
    tar[metric] = (tar[att]-1)*100
  }
  plotPerformanceSpace(P, sim, metric=metric,climData = tar,axesPercentLabel='percentage.change')
}


