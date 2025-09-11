rm(list=ls())

library(foreSIGHT)


###############
# load dates, precip, PET and streamflow data for Scott Creek
data('data_A5030502_1976_1985')
# create reference data - won't include PET here since not modelled 
clim_ref = list(times = data$times,
                P = data$P)  

###############
# create exposure space for baseline claimte - no perturbations

attPerturbType = "regGrid"
attPerturb = c('P_day_all_tot')
attPerturbSamp = c(1)
attPerturbMin = c(1)
attPerturbMax = c(1)
# will hold a number of attributes to historical values
# note: normP = P99/avg rainfall.  
attHold = c('P_day_all_P99','P_day_all_avgDSD','P_day_all_nWet',
            'P_day_DJF_tot','P_day_DJF_normP99','P_day_DJF_avgDSD','P_day_DJF_nWet',
            'P_day_MAM_tot','P_day_MAM_normP99','P_day_MAM_avgDSD','P_day_MAM_nWet',
            'P_day_JJA_tot','P_day_JJA_normP99','P_day_JJA_avgDSD','P_day_JJA_nWet',
            'P_day_SON_tot','P_day_SON_normP99','P_day_SON_avgDSD','P_day_SON_nWet')
expSpace = createExpSpace(attPerturb = attPerturb,
                               attPerturbSamp = attPerturbSamp,
                               attPerturbMin = attPerturbMin,
                               attPerturbMax = attPerturbMax,
                               attPerturbType = attPerturbType,
                               attHold = attHold)

###############
# setup model settings

modelSelection = list()

modelSelection$modelType = list()
modelSelection$modelType$P = "latent" # latent variable model

modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "seas" # parameters vary with season

modelSelection[["optimisationArguments"]] = list()
modelSelection[["optimisationArguments"]][["OFtol"]] = 0.1 # stop optimization when OF < OFtol

# set penalty weights. More weights to total, and attributes for 'all' data)
modelSelection[["penaltyAttributes"]]=c('P_day_all_tot','P_day_all_P99',
                                        'P_day_all_avgDSD','P_day_all_nWet')
modelSelection[["penaltyWeights"]] = c(3,2,2,2)

# write to JSON file
modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")
write(modelSelectionJSON, file = controlFile)

###############
# generate baseline climate scenarios

time.1 = Sys.time()
sim.base = generateScenarios(reference = clim_ref,
                                       expSpace = expSpace,
                                       controlFile = controlFile,
                                       seedID = 1,
                                       numReplicates = 1)
time.2 = Sys.time()
print(time.2-time.1)

###############
# shuffle baseline climate to introduce temporal structure at annual scale

sim.shuffle = shuffle_sim(sim=sim.base,
                          clim=clim_ref,
                          attPerturb = 'P_day_all_tot_dwellTime',
                          targetVals = c(1,1.5,2,2.5),
                          targetType = 'frac')

###############
# evaluate how changes in 'P_day_all_tot_dwellTime' affect other attributes

attSel = colnames(sim.base$expSpace$targetMat)
# other attributes
attSel = c(attSel,'P_day_all_tot_dwellTime','P_day_all_avgWSD',
           'P_day_all_P99.9','P_day_JJA_P99.9','P_day_SON_P99.9',
           'P_day_DJF_P99.9','P_day_MAM_P99.9')

att = 'P_day_all_tot_dwellTime'
par(mfrow=c(5,3),mar=c(4,7,2,1))
# plot changes in a single attribute with respect to perturbed attributes
plotPerformanceAttributesOAT(clim=clim_ref,
                             sim=sim.shuffle,
                             attPerturb=att,
                             attEval=attSel,
                             cex.main = 1.5,cex.xaxis = 1,cex.yaxis = 1)  

