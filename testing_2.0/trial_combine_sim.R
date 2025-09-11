rm(list=ls())

devtools::load_all()

clim = convert_climYMD_POSIXct(tank_obs)

############################################################

attPerturb<-c("P_day_all_tot","Temp_day_all_avg")
attPerturbType = "regGrid"
attPerturbSamp = c(2, 2)
attPerturbMin = c(0.8, -1)
attPerturbMax = c(1.1, 1)
expSpace <- createExpSpace(attPerturb = attPerturb,
                           attPerturbSamp = attPerturbSamp,
                           attPerturbMin = attPerturbMin,
                           attPerturbMax = attPerturbMax,
                           attPerturbType = attPerturbType)
data(tankDat)
simScaling <- generateScenarios(reference = clim,
                                expSpace = expSpace,
                                controlFile = "scaling")

############################################################

attPerturb<-c("P_day_all_tot")
attPerturbType = "regGrid"
attPerturbSamp = c(2)
attPerturbMin = c(0.8)
attPerturbMax = c(1.1)
expSpace.P <- createExpSpace(attPerturb = attPerturb,
                           attPerturbSamp = attPerturbSamp,
                           attPerturbMin = attPerturbMin,
                           attPerturbMax = attPerturbMax,
                           attPerturbType = attPerturbType)
simScaling.P <- generateScenarios(reference = clim,
                                expSpace = expSpace.P,
                                controlFile = "scaling")

############################################################

attPerturb<-c("Temp_day_all_avg")
attPerturbType = "regGrid"
attPerturbSamp = c(2)
attPerturbMin = c(-1)
attPerturbMax = c(1)
expSpace.Temp <- createExpSpace(attPerturb = attPerturb,
                           attPerturbSamp = attPerturbSamp,
                           attPerturbMin = attPerturbMin,
                           attPerturbMax = attPerturbMax,
                           attPerturbType = attPerturbType)
simScaling.Temp <- generateScenarios(reference = clim,
                                expSpace = expSpace.Temp,
                                controlFile = "scaling")

############################################################


sim.new = combine_sims(sim.1 = simScaling.P, var.1 = 'P',sim.2 = simScaling.Temp, var.2 = 'Temp')

