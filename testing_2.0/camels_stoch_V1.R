
rm(list=ls())

#foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT/'
foreSIGHTDir = '/scratchdata1/users/a1065639/DEW_foreSIGHT/foreSIGHT/'

devtools::load_all(foreSIGHTDir)
#devtools::install(foreSIGHTDir)
#library(foreSIGHT)

runDirname = paste0(foreSIGHTDir,'testing_2.0/')
setwd(runDirname)

############################################################################

load('data_A5050517_1976_2005.RData')
#numReplicates = 50
#cores = 50

#numReplicates = 5
numReplicates = 1
#numReplicates = 25
#cores = 25
cores = 1

#cores = parallel::detectCores()-1
#print(cores)

#clim_ref = convert_climYMD_POSIXct(barossa_obs)
#clim_ref = barossa_obs
#clim_ref$P = clim_ref$P[,1:3]
#numReplicates = 5
#numReplicates = 1
#cores = 1

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
modelSelection[["optimisationArguments"]][["nMultiStart"]] <- 2
#modelSelection[["optimisationArguments"]][["OFtol"]] <- 0.05
modelSelection[["optimisationArguments"]][['RGN.control']] = list(iterMax=5)

modelSelection[["penaltyAttributes"]] <- c('P_day_all_tot','P_day_all_avgDSD','P_day_all_P99','P_day_all_nWet','P_day_all_tot_cv')
modelSelection[["penaltyWeights"]] = rep(3,length(modelSelection[["penaltyAttributes"]]))

modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")
write(modelSelectionJSON, file = controlFile)

############################################################################

#attPerturbType = "regGrid"
#attPerturb = c('P_day_all_avgDSD')
#attPerturbSamp = c(5)
#attPerturbMin = c(1)
#attPerturbMax = c(1.4)
## attPerturbSamp = c(1)
## attPerturbMin = c(1)
## attPerturbMax = c(1)
#attHold = c('P_day_all_tot','P_day_all_P99','P_day_all_nWet','P_day_all_tot_cv')

#attPerturbType = "regGrid"
#attPerturb = c('P_day_all_tot')
#attPerturbSamp = c(5)
#attPerturbMin = c(0.8)
#attPerturbMax = c(1.2)
#attHold = c('P_day_all_avgDSD','P_day_all_P99','P_day_all_nWet','P_day_all_tot_cv')

attPerturbType = "regGrid"
attPerturb = c('P_day_all_P99')
#attPerturbSamp = c(5)
##attPerturbSamp = c(2)
#attPerturbMin = c(0.9)
#attPerturbMax = c(1.3)
attPerturbSamp = c(1)
attPerturbMin = c(1.3)
attPerturbMax = c(1.3)
##attPerturbSamp = c(1)
##attPerturbMin = c(0.9)
##attPerturbMax = c(0.9)
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

time.1 = Sys.time()
sim = generateScenarios(reference = clim_ref,
                        expSpace = expSpace,
                        controlFile = controlFile,
                        seedID = 1,
                        numReplicates = numReplicates,
                        cores = cores)
time.2 = Sys.time()
print(time.2-time.1)

##########################################################################

pdf(file=paste0(runDirname,'summary_',attPerturb,'.pdf'))

plotScenarios(sim)

#$###########################################################################
#$#
#$#PETclim = calc_ClimDaily_dayOfYearWindow(obs=clim_ref$PET,
#$#                                         dateObs = clim_ref$times,
#$#                                         dateClim = clim_ref$times,inc=14)
#$#PETclim = apply(PETclim,1,mean,na.rm=T)
#$#clim_ref_PET = clim_ref; clim_ref_PET$PET = PETclim
#$#
#$#sim.addPET = add_obs_var_to_sim(sim,var='PET',data = PETclim)

############################################################################

nTar = length(sim$Rep1)
if (nTar>1){
  
  attSel = colnames(expSpace$targetMat)
  attSel = c(attSel,'P_day_all_tot_dwellTime','P_day_all_avgWSD',
             'P_day_all_P99.9','P_day_JJA_P99.9','P_day_SON_P99.9',
             'P_day_DJF_P99.9','P_day_MAM_P99.9')
  
  Perf = calcPerformanceAttributes(clim=clim_ref,sim=sim,attSel=attSel)

print(Perf)

  par(mfrow=c(5,5),mar=c(4,5,2,1))
  for (att in names(Perf)){
    plotPerformanceOAT(Perf, sim, metric=att,plotType='base',attSel=attPerturb)
  } 
  
}

############################################################################

source(paste0(runDirname,'GR4J_funcs.R'))
source(paste0(runDirname,'boxplot.ext_DM.r'))
#devtools::load_all(foreSIGHTDir)

dates = as.Date(clim_ref$times)
Param = setup_cal_GR4J(dates = dates,P=clim_ref$P,PET=clim_ref$PET,Qobs=Qobs)
systemArgs = list(dates=dates,Param=Param,PET=clim_ref$PET)
metrics = c('meanQ','P99','P25','min3yr')

print('calc sysOutSim')
sysOutSim = runSystemModel(sim=sim,systemModel = GR4J_wrapper,systemArgs = systemArgs,metrics = metrics)
print('done') 

print('calc sysOutClim')
sysOutClim = GR4J_wrapper(data = clim_ref, systemArgs = systemArgs,metrics=metrics)
print('done')

print('calc evaluate_system_metrics')
eval = evaluate_system_metrics(sim=sim,clim=clim_ref,
                        systemModel=GR4J_wrapper,systemArgs=systemArgs,
                        metrics=metrics)
print('done')

par(mfrow=c(2,2))
for (metric in metrics){
  yAll = c(eval$systemPerf_base[[metric]],eval$systemPerf_obsClim[metric])
  ylim = c(min(yAll),max(yAll))
  boxplot.ext(as.vector(eval$systemPerf_base[[metric]]),ylim=ylim)
  points(eval$systemPerf_obsClim[metric],col='red')
  title(metric)
}

plotPerformanceOAT(performance = sysOutSim, sim=sim, metric = 'meanQ',attSel=attPerturb)
plotPerformanceOAT(performance = sysOutSim, sim=sim, metric = 'P99',attSel=attPerturb)
plotPerformanceOAT(performance = sysOutSim, sim=sim, metric = 'P25',attSel=attPerturb)
plotPerformanceOAT(performance = sysOutSim, sim=sim, metric = 'min3yr',attSel=attPerturb)

############################################################################

dev.off()

save.image(file=paste0(runDirname,'summary_',attPerturb,'.RData'))

#$#P=sim$Rep1$Target1$P$sim
#$#clim_sim = clim_ref; clim_sim$P = P
#$#
#$#print(calculateAttributes(clim_sim,'P_day_DJF_P99.9') / calculateAttributes(clim_sim,'P_day_DJF_P99'))
#$#print(calculateAttributes(clim_ref,'P_day_DJF_P99.9') / calculateAttributes(clim_ref,'P_day_DJF_P99'))
#$#
#$#
#$#P=sim$Rep1$Target2$P$sim
#$#clim_sim = clim_ref; clim_sim$P = P
#$#print(calculateAttributes(clim_sim,'P_day_DJF_P99.9') / calculateAttributes(clim_sim,'P_day_DJF_P99'))
#$#
#$#
