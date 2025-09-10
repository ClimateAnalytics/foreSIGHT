rm(list=ls())

library(foreSIGHT)

# load dates, precip, PET and streamflow data for Scott Creek
data('data_A5030502_1976_1985')

############################################################################

# create reference data - won't include PET here since not modelled 
clim_ref = list(times = data$times,
                P = data$P)  

############################################################################
# create exposure space

attPerturbType = "regGrid"
# perturb seasonality and 99% rainfall
attPerturb = c('P_day_all_seasRatioMarMay','P_day_all_P99')
# consider a large 5x5 grid
attPerturbSamp = c(5,5)
attPerturbMin = c(0.7,1.)
attPerturbMax = c(1.3,1.3)
# will hold a number of attributes to historical values
attHold = c('P_day_all_tot_m',
            'P_day_DJF_tot_m','P_day_JJA_tot_m','P_day_SON_tot_m',
            'P_day_all_nWet_m',
            'P_day_DJF_nWet_m','P_day_MAM_nWet_m',
            'P_day_JJA_nWet_m','P_day_SON_nWet_m',
            'P_day_all_avgDSD',
            'P_day_DJF_avgDSD','P_day_MAM_avgDSD',
            'P_day_JJA_avgDSD','P_day_SON_avgDSD'
            )
# tied attributes
# note: normP = P99/avg rainfall. normP in each season is tied to P99 
attTied = list(P_day_all_P99=c('P_day_DJF_normP99', 
                               'P_day_MAM_normP99',
                               'P_day_JJA_normP99',
                               'P_day_SON_normP99'),
              P_day_all_seasRatioMarMay=c('P_day_MAM_tot_m')) # MAM rain tied to seasRatioJuneAug

expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType,
                          attHold = attHold,
                          attTied = attTied)

############################################################################
# setup model settings

modelSelection = list()

modelSelection$modelType = list()
modelSelection$modelType$P = "latent" # latent variable model

modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "seas" # parameters vary with season

modelSelection[["optimisationArguments"]] = list()
modelSelection[["optimisationArguments"]][["OFtol"]] = 0.1 # stop optimization when OF < OFtol

# set penalty weights. More weights to perturbed attributes and total rainfall (which is known to have large impact)
modelSelection[["penaltyAttributes"]]=c('P_day_all_seasRatioMarMay',
                                        'P_day_all_P99','P_day_all_tot_m',
                                        'P_day_all_avgDSD','P_day_all_nWet_m')
modelSelection[["penaltyWeights"]] = c(3,3,3,1.5,1.5)

# write to JSON file
modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")
write(modelSelectionJSON, file = controlFile)

###############
# generate scenarios

time.1 = Sys.time()
sim.stoch = generateScenarios(reference = clim_ref,
                                       expSpace = expSpace,
                                       controlFile = controlFile,
                                       seedID = 1,
                                       numReplicates = 50,
                                       cores = 25)
time.2 = Sys.time()
print(time.2-time.1)

save(file='../egScottCreekSimStoch.rda',sim.stoch)
