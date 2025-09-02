rm(list=ls())

library(foreSIGHT)

#devtools::load_all()

timeStart = as.POSIXct('2000/01/01 00:00:00',tz='UTC')
timeEnd = as.POSIXct('2000/12/31 23:00:00',tz='UTC')

# timeStart = as.POSIXct('2000/01/01 00:00:00',tz='UTC')
# timeEnd = as.POSIXct('2009/12/31 23:00:00',tz='UTC')

times = seq(timeStart,timeEnd,by='hours')

nTimes = length(times)

lambda = 0.02
gamma = 1/10
beta = 0.3
eta = 2
mux = 4
t.sim = nTimes
seed = 1 
set.seed(seed)

simulation = SWGsim.BLRPM(SWGpar=list(lambda=lambda,gamma=gamma,beta=beta,eta=eta,mux=mux),
                           nTimes=nTimes,
                           randomTerm=list(seed=seed))

subdaily_synthetic_obs = list(times=times,
                              P=simulation)

save(file='data/subdailySyntheticDat.rda',subdaily_synthetic_obs)

########################################

modelSelection = list()
modelSelection$modelType = list()
modelSelection$modelType$P = "BLRPM"
modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "ann"

modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)

controlFile = paste0(tempdir(), "\\eg_controlFile.json")
write(modelSelectionJSON, file = controlFile)

########################################

attPerturb = c('P_day_all_tot')
attHold = c('P_hour_all_sd','P_hour_all_cor','P_hour_all_nWet',
            'P_3hour_all_sd','P_3hour_all_cor','P_3hour_all_nWet',
            'P_12hour_all_sd','P_12hour_all_cor','P_12hour_all_nWet',
            'P_day_all_sd','P_day_all_cor','P_day_all_nWet')

attPerturbType = "regGrid"
attPerturbSamp = c(1)
attPerturbMin = c(1.)
attPerturbMax = c(1.)

# create the exposure space
expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType,
                          attHold = attHold)

########################################

sim = generateScenarios(reference = subdaily_synthetic_obs,
                        expSpace = expSpace,
                        controlFile = controlFile,
                        seedID = 1)

plotScenarios(sim)


