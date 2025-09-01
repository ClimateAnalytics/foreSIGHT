rm(list=ls())

foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT/'
devtools::load_all(foreSIGHTDir)

#library(foreSIGHT)

# Example 5: Multi-site stochastic simulation
#-----------------------------------------------------------------------
attPerturb <- c("P_day_all_tot_m")
attHold <- c(
  "P_day_all_wettest6monSeasRatio", "P_day_all_wettest6monPeakDay",
  "P_day_all_P99", "P_day_all_avgWSD_m", "P_day_all_nWetT0.999_m"
)
attPerturbType <- "regGrid"
# consider unperturbed climates in this example
attPerturbSamp <- attPerturbMin <- attPerturbMax <- c(1)
expSpace <- createExpSpace(
  attPerturb = attPerturb,
  attPerturbSamp = attPerturbSamp,
  attPerturbMin = attPerturbMin,
  attPerturbMax = attPerturbMax,
  attPerturbType = attPerturbType,
  attHold = attHold
)
# load multi-site rainfall data
data(barossaDat)
clim_ref = barossa_obs; clim_ref$P = clim_ref$P[,1:3]
# specify the penalty settings in a list
controlFileList <- list()
controlFileList[["penaltyAttributes"]] <- c(
  "P_day_all_tot_m",
  "P_day_all_wettest6monSeasRatio", "P_day_all_wettest6monPeakDay"
)
controlFileList[["penaltyWeights"]] <- c(0.5, 0.5, 0.5)
# specify the alternate model selections
controlFileList[["modelType"]] <- list()
controlFileList[["modelType"]][["P"]] <- "latent"
# specify model parameter selection
controlFileList[["modelParameterVariation"]] <- list()
controlFileList[["modelParameterVariation"]][["P"]] <- "har"
# specify settings for multi-site model
controlFileList[["spatialOptions"]] <- list()
# specify spatial correlation perturbation factor
#controlFileList[["spatialOptions"]][["spatCorFac"]] <- 0.9
# specify optimization arguments
controlFileList[["optimisationArguments"]] <- list()
controlFileList[["optimisationArguments"]][["OFtol"]] <- 0.02
# write control file settings to file
controlFileJSON <- jsonlite::toJSON(controlFileList, pretty = TRUE, auto_unbox = TRUE)
write(controlFileJSON, file = paste0(tempdir(), "controlFile.json"))
# run multi-site stochastic simulation - this will take a long time (e.g. hours)
egMultiSiteSim <- generateScenarios(
  reference = clim_ref, expSpace = expSpace,
  controlFile = paste0(tempdir(), "controlFile.json"), seed = 1
)
save(file='data/egMultiSiteSim.rda',egMultiSiteSim)


