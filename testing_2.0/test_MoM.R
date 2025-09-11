rm(list=ls())

devtools::load_all()

data("tankDat")

calculateAttributes(tank_obs,attSel='P_day_all_tot_m')

#modelTag = 'P-ann-wgen'
modelTag = 'P-seas-wgen'

par = modCalibrator(obs=tank_obs,modelTag = modelTag)

dates = tank_obs$times
# datStart = dates[1]
# datFinish = dates[length(dates)]

#sim = modSimulator(datStart=datStart,datFinish=datFinish,parS = par,modelTag = modelTag)

datInd = list(day=get.date.ind(times=dates))

seed=1
set.seed(seed)
randomVector <- stats::runif(n=length(dates)) # Random vector to be passed into weather generator to reduce runtime


P = simClim(parS=par[[modelTag]],modelTag=modelTag,
        modelInfo=modelInfoList[[modelTag]],
        datInd=datInd,
        randomTerm=list(randomVector=randomVector,seed=seed))


clim_sim = tank_obs
clim_sim$P = P
par1 = modCalibrator(obs=clim_sim,modelTag = modelTag)



data(tankDat) #Load tank data (tank_obs)
modelTag=c("P-ann-wgen","Temp-har-wgenLM") #Select a rainfall and a temperature generator
out<- modCalibrator(obs = tank_obs, #Calibrate models
                    modelTag = modelTag)

