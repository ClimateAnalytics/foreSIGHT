rm(list=ls())

devtools::load_all()

# # Example 1: Simple scaling
# #-----------------------------------------------------------------------
# attPerturb<-c("P_ann_tot_m","Temp_ann_avg_m")
# # attPerturb<-c("P_day_all_tot_m","Temp_day_all_avg_m") # this could also be "P_year_all_avg"
# # attPerturb<-c("P_year_all_avg","Temp_day_all_avg")
# attPerturbType = "regGrid"
# attPerturbSamp = c(2, 2)
# attPerturbMin = c(0.8, -1)
# attPerturbMax = c(1.1, 1)
# expSpace <- createExpSpace(attPerturb = attPerturb,
#                            attPerturbSamp = attPerturbSamp,
#                            attPerturbMin = attPerturbMin,
#                            attPerturbMax = attPerturbMax,
#                            attPerturbType = attPerturbType)
# data(tankDat)
# simScaling <- generateScenarios(reference = tank_obs,
#                                 expSpace = expSpace,
#                                 controlFile = "scaling")
# 
# 
# # checks
# Psim = simScaling$Rep1$Target1$P
# mean(Psim)/mean(tank_obs$P) # should be 0.8
# 
# Tsim = simScaling$Rep1$Target1$Temp
# mean(Tsim) - mean(tank_obs$Temp) # should be -1

# Example 2: Seasonal scaling
#-----------------------------------------------------------------------
# attPerturb<-c("P_ann_tot_m","P_ann_seasRatio")
# # attPerturb<-c("P_day_all_tot_m","P_day_all_seasRatio")
# #attPerturb<-c("P_year_all_avg","P_day_all_seasRatio")
# attPerturbType = "regGrid"
# attPerturbSamp = c(2, 2)
# attPerturbMin = c(0.8, 0.9)
# attPerturbMax = c(1.1, 1.2)
# expSpace <- createExpSpace(attPerturb = attPerturb,
#                            attPerturbSamp = attPerturbSamp,
#                            attPerturbMin = attPerturbMin,
#                            attPerturbMax = attPerturbMax,
#                            attPerturbType = attPerturbType)
# data(tankDat)
# seasScaling <- generateScenarios(reference = tank_obs,
#                                  expSpace = expSpace,
#                                  controlFile = "scaling")
# 
# # checks
# clim_sim = tank_obs; clim_sim$P = seasScaling$Rep1$Target1$P
# calculateAttributes(clim_sim,attSel=attPerturb)/calculateAttributes(tank_obs,attSel=attPerturb)
# 
# # Psim = seasScaling$Rep1$Target1$P
# # mean(Psim)/mean(tank_obs$P)
# 
# # #keep = which(as.integer(format(tank_obs$times,'%m'))%in%(seq(3,8)))
# # keep = which(as.integer(format(tank_obs$times,'%m'))%in%(seq(3,8)))
# # SRsim = mean(Psim[keep])/mean(Psim)
# # SRobs = mean(tank_obs$P[keep])/mean(tank_obs$P)
# # SRsim/SRobs
# # 
# # pause

# Example 3: Stochastic simulation using foreSIGHT default settings
#----------------------------------------------------------------------
# create an exposure space
# attPerturb <- c("P_ann_tot_m", "P_ann_nWet_m", "P_ann_R10_m")
# attHold <- c("P_Feb_tot_m", "P_SON_dyWet_m", "P_JJA_avgWSD_m", "P_MAM_tot_m",
#             "P_DJF_avgDSD_m", "Temp_ann_rng_m", "Temp_ann_avg_m")
# # attPerturb <- c("P_day_all_tot_m", "P_day_all_nWet_m", "P_day_all_R10_m")
# # attHold <- c("P_day_Feb_tot_m", "P_day_SON_dyWet_m", "P_day_JJA_avgWSD_m", "P_day_MAM_tot_m",
# #              "P_day_DJF_avgDSD_m", "Temp_day_all_rng_m", "Temp_day_all_avg_m")
# attPerturbType = "regGrid"
# attPerturbSamp = c(2, 1, 1)
# attPerturbMin = c(0.8, 1, 1)
# attPerturbMax = c(1.1, 1, 1)
# 
# # attPerturb <- c("Temp_day_all_rng_m")
# # attHold <- c("Temp_day_all_avg_m")
# # attPerturbType = "regGrid"
# # attPerturbSamp = c(2)
# # attPerturbMin = c(0.8)
# # attPerturbMax = c(1.1)
# 
# expSpace <- createExpSpace(attPerturb = attPerturb,
#                            attPerturbSamp = attPerturbSamp,
#                            attPerturbMin = attPerturbMin,
#                            attPerturbMax = attPerturbMax,
#                            attPerturbType = attPerturbType,
#                            attHold = attHold)
# # load example data available in foreSIGHT
# data(tankDat)
# # perform stochastic simulation
# simStochastic <- generateScenarios(reference = tank_obs,
#                                    expSpace = expSpace,
#                                    simLengthNyrs = 30,seedID=1)
# 
# plotScenarios(simStochastic)
# 
# pause

# ## End(Not run)

# Example 4: Simple Scaling with multi-site data
#-----------------------------------------------------------------------
# attPerturb <- c("P_ann_tot_m","P_ann_seasRatio")
# # attPerturb <- c("P_day_all_tot_m","P_day_all_seasRatio")
# attPerturbType = "regGrid"
# attPerturbSamp = c(3, 3)
# attPerturbMin = c(0.8, 0.8)
# attPerturbMax = c(1.2, 1.2)
# expSpace <- createExpSpace(attPerturb = attPerturb,
#                            attPerturbSamp = attPerturbSamp,
#                            attPerturbMin = attPerturbMin,
#                            attPerturbMax = attPerturbMax,
#                            attPerturbType = attPerturbType)
# # load multi-site rainfall data
# data(barossaDat)
# # perform simple scaling
# simScaling <- generateScenarios(reference = barossa_obs,
#                                 expSpace = expSpace,
#                                 controlFile = "scaling")
# # checks
# # mean rain at site=1 and target=1 - expect 0.8
# mean(simScaling$Rep1$Target1$P[,1]) / mean(barossa_obs$P[,1])
# # mean rain at site=10 and target=2 - expect 1.2
# mean(simScaling$Rep1$Target3$P[,1]) / mean(barossa_obs$P[,1])
# 
# pause

# Example 5: Multi-site stochastic simulation
#-----------------------------------------------------------------------
## Not run: 
attPerturb <- c("P_ann_tot_m")
attHold <- c("P_ann_wettest6monSeasRatio","P_ann_wettest6monPeakDay",
            "P_ann_P99","P_ann_avgWSD_m", "P_ann_nWetT0.999_m")
# attPerturb <- c("P_day_all_tot_m")
# attHold <- c("P_day_all_wettest6monSeasRatio","P_day_all_wettest6monPeakDay",
#              "P_day_all_P99","P_day_all_avgWSD_m", "P_day_all_nWetT0.999_m")
attPerturbType = "regGrid"
# consider unperturbed climates in this example
attPerturbSamp = attPerturbMin = attPerturbMax = c(1)
expSpace <- createExpSpace(attPerturb = attPerturb,
                           attPerturbSamp = attPerturbSamp,
                           attPerturbMin = attPerturbMin,
                           attPerturbMax = attPerturbMax,
                           attPerturbType = attPerturbType,
                           attHold = attHold)
# load multi-site rainfall data
data(barossaDat)

clim = barossa_obs
# clim$P = clim$P[,1:2]

# specify the penalty settings in a list
controlFileList <- list()
controlFileList[["penaltyAttributes"]] <- c("P_ann_tot_m",
                                           "P_ann_wettest6monSeasRatio","P_ann_wettest6monPeakDay")
# controlFileList[["penaltyAttributes"]] <- c("P_day_all_tot_m",
#                                             "P_day_all_wettest6monSeasRatio","P_day_all_wettest6monPeakDay")
controlFileList[["penaltyWeights"]] <- c(0.5,0.5,0.5)
# specify the alternate model selections
controlFileList[["modelType"]] <- list()
controlFileList[["modelType"]][["P"]] <- "latent"
# controlFileList[["modelType"]][["P"]] <- "LV"
# specify model parameter selection
controlFileList[["modelParameterVariation"]] <- list()
controlFileList[["modelParameterVariation"]][["P"]] <- "harmonic"
# controlFileList[["modelParameterVariation"]][["P"]] <- "har"
# specify settings for multi-site model
controlFileList[["spatialOptions"]] <- list()
# specify spatial correlation perturbation factor
controlFileList[["spatialOptions"]][["spatCorFac"]] = 0.9
# write control file sttings to file
controlFileJSON <- jsonlite::toJSON(controlFileList, pretty = TRUE, auto_unbox = TRUE)
write(controlFileJSON, file = paste0(tempdir(), "controlFile.json"))
# run multi-site stochastic simulation - this will take a long time (e.g. hours)
sim <- generateScenarios(reference = clim, expSpace = expSpace,
                         controlFile = paste0(tempdir(), "controlFile.json"),seed=1)
## End(Not run)

