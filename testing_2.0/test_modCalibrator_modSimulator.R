rm(list=ls())

devtools::load_all()

data(tankDat); obs=tank_obs                     #Get observed data

modelTagList = list()

modelTagList[[1]] = "P-ann-wgen"
modelTagList[[2]] = "P-seas-wgen"
modelTagList[[3]] = "P-har-wgen"
modelTagList[[4]] = "Temp-har-wgenO"
modelTagList[[5]] = c("P-seas-wgen","Temp-harWD-wgenO")

# modelTagList[[1]] = c("P-seas-wgen","Temp-harWD-wgenO")

for (m in 1:length(modelTagList)){
  modelTag = modelTagList[[m]]
  pars=modCalibrator(obs=obs,modelTag=modelTag)   #Calibrate models
  sim=modSimulator(datStart="2007-01-01",         #Simulate!
                   datFinish="2106-12-31",
                   modelTag=modelTag,
                   parS=pars,
                   seed=123,
                   file=paste0("tester.csv"),
                   IOmode="verbose")
  pars.1=modCalibrator(obs=sim,modelTag=modelTag)   #Calibrate models
  print(modelTag)
  print(pars)
  print(pars.1)
}