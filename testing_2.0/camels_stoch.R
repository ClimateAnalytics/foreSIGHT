# # scaling
# 
# settings.R
# 
# load data camels
# 
# expSpace
# 
# sim scaling
# 
# system response
# 
# plot OAT, 2d
# 
# 
# # stoch
# 
# settings.R
# 
# load data camels
# 
# expSpace - OAT (Ptot, P99, DSD) (hold - nwet, seas, ann CV)
# 
# model settings
# 
# sim stoch P (not shuffle)
# sim stoch P shuffle
# sim scaling PET
# 
# combine
# 
# system response
# 
# plot OAT


rm(list=ls())

source(paste0(runDirname,'settings.R'))

runDirname = paste0(foreSIGHTDir,'testing_2.0/')

############################################################################
# read catchment data

#catchment = 'A5030502' # Scott Creek
# catchment = 'A5050517' # North Para River at Penrice (A5050517)
# 
# startYr = 1976
# endYr = 2005

# data = load_camels(catchment,startYr,endYr)
# Qobs = data$Qobs
# 
# # create reference data 
# clim_ref <- list(times = data$times,
#                  P = data$P,
#                  PET = data$PET)

#fname = paste0(runDirname,'data_',catchment,'_',startYr,'_',endYr,'.RData')
# save(file=fname,clim_ref,Qobs)

load('C:/Users/a1065639/Work/foreSIGHT/testing_2.0/data_A5050517_1976_2005.RData')

############################################################################
# setup model settings

modelSelection = list()

modelSelection$modelType = list()
modelSelection$modelType$P = "latent"

modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "seas"

modelSelection$postProcessing = list(P=list())
modelSelection$postProcessing$P$types = c('annVar','scaleExtremesSeas')

modelSelection[["optimisationArguments"]] <- list()
modelSelection[["optimisationArguments"]][["nMultiStart"]] <- 5
modelSelection[["optimisationArguments"]][["OFtol"]] <- 0.05

modelSelection[["penaltyAttributes"]] <- c('P_day_all_tot','P_day_all_avgDSD','P_day_all_P99','P_day_all_nWet','P_day_all_tot_cv')
modelSelection[["penaltyWeights"]] = rep(3,length(modelSelection[["penaltyAttributes"]]))

modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")
write(modelSelectionJSON, file = controlFile)

# ############################################################################
# # setup baseline calibration - used for faster parameter estimation for perturbed climates
# 
# attPerturbType = "regGrid"
# attPerturb = c('P_day_all_avgDSD')
# attPerturbSamp = c(1)
# attPerturbMin = c(1)
# attPerturbMax = c(1)
# attHold = c('P_day_all_tot','P_day_all_P99','P_day_all_nWet','P_day_all_tot_cv')
# attsAll = c(attPerturb,attHold)
# 
# expSpace = createExpSpace(attPerturb = attPerturb,
#                           attPerturbSamp = attPerturbSamp,
#                           attPerturbMin = attPerturbMin,
#                           attPerturbMax = attPerturbMax,
#                           attPerturbType = attPerturbType,
#                           attHold = attHold)
# 
# attTied = list(seas=attsAll[!attsAll%in%c('P_day_all_tot_cv')])
# expSpace = tieAttributes(expSpace=expSpace,attTied=attTied)
# 
# ############################################################################
# 
# time.1 = Sys.time()
# sim = generateScenarios(reference = clim_ref,
#                               expSpace = expSpace,
#                               controlFile = controlFile,
#                               seedID = 1,
#                               numReplicates = 5,
#                               cores = 1)
# time.2 = Sys.time()
# print(time.2-time.1)
# 
# plotScenarios(sim)
# 
# ############################################################################
# 
# modelSelection[["optimisationArguments"]]$suggestions = sim$Rep1$Target1$P$par
# modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
# controlFile = paste0(tempdir(), "\\eg_controlFile.json")
# write(modelSelectionJSON, file = controlFile)
# 

############################################################################

attPerturbType = "regGrid"
attPerturb = c('P_day_all_avgDSD')
attPerturbSamp = c(5)
attPerturbMin = c(1)
attPerturbMax = c(1.4)
# attPerturbSamp = c(1)
# attPerturbMin = c(1)
# attPerturbMax = c(1)
attHold = c('P_day_all_tot','P_day_all_P99','P_day_all_nWet','P_day_all_tot_cv')
#attHold = c('P_day_all_tot','P_day_all_xP99overPave','P_day_all_nWet','P_day_all_tot_cv')


# attPerturbType = "regGrid"
# attPerturb = c('P_day_all_tot')
# attPerturbSamp = c(5)
# attPerturbMin = c(0.8)
# attPerturbMax = c(1.2)
# # attPerturbSamp = c(1)
# # attPerturbMin = c(1)
# # attPerturbMax = c(1)
# attHold = c('P_day_all_avgDSD','P_day_all_P99','P_day_all_nWet','P_day_all_tot_cv')
# #attHold = c('P_day_all_tot','P_day_all_xP99overPave','P_day_all_nWet','P_day_all_tot_cv')

attPerturbType = "regGrid"
attPerturb = c('P_day_all_P99')
attPerturbSamp = c(5)
attPerturbMin = c(0.9)
attPerturbMax = c(1.3)
attHold = c('P_day_all_tot','P_day_all_avgDSD','P_day_all_nWet','P_day_all_tot_cv')

#############################

attsAll = c(attPerturb,attHold)

expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType,
                          attHold = attHold)

attTied = list(seas=attsAll[!attsAll%in%c('P_day_all_tot_cv')])
expSpace = tieAttributes(expSpace=expSpace,attTied=attTied)

############################################################################

# time.1 = Sys.time()
# sim.new = generateScenarios(reference = clim_ref,
#                         expSpace = expSpace,
#                         controlFile = controlFile,
#                         seedID = 1,
#                         numReplicates = numReplicates,
#                         cores = cores)
# time.2 = Sys.time()
# print(time.2-time.1)
# 
# plotScenarios(sim.new)
# 
# ##########################################################################
# 
# sim.new.shuffle = shuffle_sim(sim=sim.new,clim = clim_ref,
#                               attPerturb = 'P_day_all_tot_dwellTime',
#                               targetVals = 0)


sim.new.shuffle = shuffle_sim(sim=sim,clim = clim_ref,
                              attPerturb = 'P_day_all_tot_dwellTime',
                              targetVals = c(1,1.5,2))

##########################################################################

PETclim = calc_ClimDaily_dayOfYearWindow(obs=clim_ref$PET,
                                         dateObs = clim_ref$times,
                                         dateClim = clim_ref$times,inc=14)

PETclim = apply(PETclim,1,mean,na.rm=T)

sim.new.shuffle.addPET = add_obs_var_to_sim(sim.new.shuffle,var='PET',data = PETclim)

############################################################################


attSel = colnames(expSpace$targetMat)
attSel = c(attSel,'P_day_all_tot_dwellTime','P_day_all_P99.9','P_day_all_avgWSD')

P = calcPerformanceAttributes(clim=clim_ref,sim=sim.new.shuffle.addPET,attSel=attSel)

attPerturb = 'P_day_all_tot_dwellTime'

par(mfrow=c(5,5),mar=c(4,5,2,1))
for (att in names(P)){
  plotPerformanceOAT(P, sim.new.shuffle.addPET, metric=att,plotType='base',attSel=attPerturb)
}

############################################################################

source('testing_2.0/GR4J_funcs.R')
  
dates = as.Date(data$times)

Param = setup_cal_GR4J(dates = dates,data$P,data$PET,data$Qobs)

systemArgs = list(dates=dates,Param=Param)

sysOutSim = runSystemModel(sim=sim.new.shuffle.addPET,systemModel = GR4J_wrapper,systemArgs = systemArgs,metrics = c('meanQ','P99','P25','min3yr'),varNames=c('P','PET'))
  
clim_ref_PET = clim_ref; clim_ref_PET$PET = PETclim

sysOutClim = GR4J_wrapper(data = clim_ref_PET, systemArgs = systemArgs,metrics = c('meanQ','P99','P25','min3yr'))
  
plotPerformanceOAT(performance = sysOutSim, sim=sim.new.shuffle.addPET, metric = 'meanQ',attSel=attPerturb)
plotPerformanceOAT(performance = sysOutSim, sim=sim.new.shuffle.addPET, metric = 'P99',attSel=attPerturb)
plotPerformanceOAT(performance = sysOutSim, sim=sim.new.shuffle.addPET, metric = 'P25',attSel=attPerturb)
plotPerformanceOAT(performance = sysOutSim, sim=sim.new.shuffle.addPET, metric = 'min3yr',attSel=attPerturb)

#plotPerformanceSpace(performance = sysOutSim, sim=sim.new.shuffle.addPET, metric = 'meanQ')


