rm(list=ls())

# setwd('C:/Users/a1065639/Work/foreSIGHT/')
# attPerturb <- c("P_day_all_tot_m","P_day_all_seasRatio")
# attHold <- c("P_day_all_nWet_m", "P_day_all_R10_m", "P_day_Feb_tot_m", "P_day_SON_dyWet_m", "P_day_JJA_avgWSD_m", "P_day_MAM_tot_m", "P_day_DJF_avgDSD_m", "Temp_day_all_rng_m", "Temp_day_all_avg_m")
# controlFile = "C:/Users/a1065639/Work/foreSIGHT/testing_2.0/controlFile_brees_paper_foreSIGHTnew.json"

setwd('C:/Users/a1065639/Work/foreSIGHT1.2/foreSIGHT/')
attPerturb <- c("P_ann_tot_m","P_ann_seasRatio")
attHold <- c("P_ann_nWet_m", "P_ann_R10_m", "P_Feb_tot_m", "P_SON_dyWet_m", "P_JJA_avgWSD_m", "P_MAM_tot_m", "P_DJF_avgDSD_m", "Temp_ann_rng_m", "Temp_ann_avg_m")
controlFile = "C:/Users/a1065639/Work/foreSIGHT/testing_2.0/controlFile_brees_paper.json.txt"

devtools::load_all()

attPerturbType <- "regGrid"
attPerturbSamp <- c(10,16)
attPerturbMin <- c(0.8, 0.80)
attPerturbMax <- c(1.1, 1.3)
expSpace <- createExpSpace(attPerturb = attPerturb, attPerturbSamp = attPerturbSamp, attPerturbMin = attPerturbMin, attPerturbMax = attPerturbMax, attPerturbType = attPerturbType, attHold = attHold)

expSpace_subset <- expSpace
expSpace_subset$targetMat <- expSpace$targetMat[1,]
data("tankDat")
simTest <- generateScenarios(reference = tank_obs, 
                             expSpace = expSpace_subset, 
                             simLengthNyrs = 300, 
                             numReplicates = 1,
                             seedID = 2,
                             controlFile = controlFile)

plotScenarios(simTest)