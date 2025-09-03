# simulations similar to those in Bree's foreSIGHT paper

rm(list=ls())

controlFile = "egSim_controlFile.json"

library(foreSIGHT)

############################################################################################################

attPerturb <- c("P_day_all_tot_m","P_day_all_seasRatio")
attHold <- c("P_day_all_nWet_m", "P_day_all_R10_m", "P_day_Feb_tot_m", "P_day_SON_dyWet_m", "P_day_JJA_avgWSD_m", "P_day_MAM_tot_m", "P_day_DJF_avgDSD_m", "Temp_day_all_rng_m", "Temp_day_all_avg_m")
attPerturbType <- "regGrid"
attPerturbSamp <- c(7,7)
attPerturbMin <- c(0.8, 0.9)
attPerturbMax <- c(1.1, 1.2)
expSpace <- createExpSpace(attPerturb = attPerturb, attPerturbSamp = attPerturbSamp, attPerturbMin = attPerturbMin, attPerturbMax = attPerturbMax, attPerturbType = attPerturbType, attHold = attHold)

data("tankDat")
time.1 = Sys.time()
simTest <- generateScenarios(reference = tank_obs, 
                             expSpace = expSpace, 
                             simLengthNyrs = 100, 
                             numReplicates = 10,
                             seedID = 1,
                             controlFile = controlFile,
			     cores=25)
time.2 = Sys.time()
print(time.2-time.1)

egSimSummary = getSimSummary(simTest)
save(file='../egSimSummary.rda',egSimSummary)

############################################################################################################

metrics = c("average daily deficit (L)", "reliability (fraction)")

systemArgs <- list(roofArea = 205,                         # roof area in m2
                    nPeople = 1,                            # number of people using water
                    tankVol = 2400,                         # tank volume in L
                    firstFlush = 2.0,                       # first flush volume
                    write.file = F,                         # write output tank timeseries to file T/F?
                    metrics = metrics
		    )

egSimPerformance <- runSystemModel(sim=simTest,
                              systemModel=tankWrapper,
                              systemArgs=systemArgs,
                              metrics = metrics)
save(file='../egSimPerformance.rda',egSimPerformance)

############################################################################################################

systemArgsB <- list(roofArea = 205,                         # roof area in m2
                    nPeople = 1,                            # number of people using water
                    tankVol = 2600,                         # tank volume in L
                    firstFlush = 2.0,                       # first flush volume
                    write.file = F,                         # write output tank timeseries to file T/F?
                    metrics = metrics
                    )

egSimPerformanceB <- runSystemModel(sim=simTest,
                              systemModel=tankWrapper,
                              systemArgs=systemArgsB,
                              metrics = metrics)
save(file='../egSimPerformanceB.rda',egSimPerformanceB)

