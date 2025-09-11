rm(list=ls())

# foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT/'
# devtools::load_all(foreSIGHTDir)

library(foreSIGHT)

attPerturb <- c("P_day_all_tot_m", "Temp_day_all_avg")
attPerturbType <- "regGrid"
attPerturbSamp <- c(10, 10)
attPerturbMin <- c(0.8, -1)
attPerturbMax <- c(1.2, 1)
expSpace <- createExpSpace(
  attPerturb = attPerturb,
  attPerturbSamp = attPerturbSamp,
  attPerturbMin = attPerturbMin,
  attPerturbMax = attPerturbMax,
  attPerturbType = attPerturbType
)
data(tankDat)
simScaling <- generateScenarios(
  reference = tank_obs,
  expSpace = expSpace,
  controlFile = "scaling"
)
egScalSummary = getSimSummary(simScaling)

# use the simulation to run a system model
systemArgs <- list(
  roofArea = 205, nPeople = 1, tankVol = 2400,
  firstFlush = 2.0, write.file = FALSE
)
tankMetrics <- viewTankMetrics()
egScalPerformance <- runSystemModel(
  sim = simScaling,
  systemModel = tankWrapper,
  systemArgs = systemArgs,
  metrics = tankMetrics[1:2]
)

save(file='data/egScalSummary.rda',egScalSummary)
save(file='data/egScalPerformance.rda',egScalPerformance)

