rm(list=ls())

setwd('C:/Users/a1065639/Work/foreSIGHT/')
controlFile = "C:/Users/a1065639/Work/foreSIGHT/testing_2.0/controlFile_brees_paper_foreSIGHTnew.json"

# setwd('C:/Users/a1065639/Work/foreSIGHT1.2/foreSIGHT/')
# attPerturb <- c("P_ann_tot_m","P_ann_seasRatio")
# attHold <- c("P_ann_nWet_m", "P_ann_R10_m", "P_Feb_tot_m", "P_SON_dyWet_m", "P_JJA_avgWSD_m", "P_MAM_tot_m", "P_DJF_avgDSD_m", "Temp_ann_rng_m", "Temp_ann_avg_m")
# controlFile = "C:/Users/a1065639/Work/foreSIGHT/testing_2.0/controlFile_brees_paper.json.txt"

devtools::load_all()

attPerturb <- c("P_day_all_tot_m","P_day_all_seasRatio")
attHold <- c("P_day_all_nWet_m", "P_day_all_R10_m", "P_day_Feb_tot_m", "P_day_SON_dyWet_m", "P_day_JJA_avgWSD_m", "P_day_MAM_tot_m", "P_day_DJF_avgDSD_m", "Temp_day_all_rng_m", "Temp_day_all_avg_m")
attPerturbType <- "regGrid"
attPerturbSamp <- c(13,13)
#attPerturbMin <- c(0.8, 0.80)
#attPerturbMax <- c(1.1, 1.3)
attPerturbMin <- c(0.8, 0.9)
attPerturbMax <- c(1.1, 1.2)
expSpace <- createExpSpace(attPerturb = attPerturb, attPerturbSamp = attPerturbSamp, attPerturbMin = attPerturbMin, attPerturbMax = attPerturbMax, attPerturbType = attPerturbType, attHold = attHold)

# attPerturb <- c("P_day_all_tot_m","P_day_all_seasRatio","P_day_all_nWet_m", "P_day_all_R10_m")
# attHold = c("P_day_Feb_tot_m", "P_day_SON_dyWet_m", "P_day_JJA_avgWSD_m", "P_day_MAM_tot_m", "P_day_DJF_avgDSD_m", "Temp_day_all_rng_m", "Temp_day_all_avg_m")
# attPerturbType <- "OAT"
# attPerturbSamp <- c(5,5,5,5)
# attPerturbMin <- c(0.8,0.8,0.85,0.9)
# attPerturbMax <- c(1.1, 1.3,1.05,1.25)
# expSpace <- createExpSpace(attPerturb = attPerturb, attPerturbSamp = attPerturbSamp, attPerturbMin = attPerturbMin, attPerturbMax = attPerturbMax, attPerturbType = attPerturbType, attHold = attHold)


expSpace_subset <- expSpace
expSpace_subset$targetMat <- expSpace$targetMat[1,]
data("tankDat")
time.1 = Sys.time()
simTest <- generateScenarios(reference = tank_obs, 
                             expSpace = expSpace_subset, 
                             simLengthNyrs = 300, 
                             numReplicates = 1,
                             seedID = 2,
                             controlFile = controlFile)
time.2 = Sys.time()
print(time.2-time.1)

plotScenarios(simTest)