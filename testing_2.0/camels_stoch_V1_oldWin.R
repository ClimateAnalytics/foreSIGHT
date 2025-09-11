
rm(list=ls())

foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT/'
devtools::load_all(foreSIGHTDir)

#library(foreSIGHT)

runDirname = paste0(foreSIGHTDir,'testing_2.0/')
setwd(runDirname)

############################################################################

load('data_A5050517_1976_2005.RData')
numReplicates = 1
cores = 1
#cores=2
#numReplicates = 5

# clim_ref = convert_climYMD_POSIXct(barossa_obs)
# clim_ref$P = clim_ref$P[,1:2]

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
modelSelection[["optimisationArguments"]][["nMultiStart"]] <- 5
modelSelection[["optimisationArguments"]][["OFtol"]] <- 0.05
#modelSelection[["optimisationArguments"]][['RGN.control']] = list(iterMax=100)
#modelSelection[["optimisationArguments"]][['RGN.control']] = list(iterMax=10)

modelSelection[["penaltyAttributes"]] <- c('P_day_all_tot','P_day_all_avgDSD','P_day_all_P99','P_day_all_nWet','P_day_all_tot_cv')
modelSelection[["penaltyWeights"]] = rep(3,length(modelSelection[["penaltyAttributes"]]))

modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")
write(modelSelectionJSON, file = controlFile)

############################################################################

# attPerturbType = "regGrid"
# attPerturb = c('P_day_all_avgDSD')
# attPerturbSamp = c(3)
# attPerturbMin = c(1)
# attPerturbMax = c(1.4)
# # attPerturbSamp = c(1)
# # attPerturbMin = c(1)
# # attPerturbMax = c(1)
# attHold = c('P_day_all_tot','P_day_all_P99','P_day_all_nWet','P_day_all_tot_cv')

# attPerturbType = "regGrid"
# attPerturb = c('P_day_all_tot')
# attPerturbSamp = c(2)
# attPerturbMin = c(0.8)
# attPerturbMax = c(1.2)
# # attPerturbSamp = c(1)
# # attPerturbMin = c(1)
# # attPerturbMax = c(1)
# attHold = c('P_day_all_avgDSD','P_day_all_P99','P_day_all_nWet','P_day_all_tot_cv')

attPerturbType = "regGrid"
attPerturb = c('P_day_all_P99')
attPerturbSamp = c(2)
attPerturbMin = c(0.9)
attPerturbMax = c(1.3)
attHold = c('P_day_all_tot','P_day_all_avgDSD','P_day_all_nWet','P_day_all_tot_cv')

#############################


expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType,
                          attHold = attHold)

attsAll = c(attPerturb,attHold)
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

#pdf(file=paste0(runDirname,'summary_',attPerturb,'.pdf'))

#plotScenarios(sim)

##########################################################################

nTar = length(sim$Rep1)
if (nTar>1){
  
  attSel = colnames(expSpace$targetMat)
  attSel = c(attSel,'P_day_all_tot_dwellTime','P_day_all_avgWSD',
             'P_day_all_P99.9','P_day_JJA_P99.9','P_day_SON_P99.9',
             'P_day_DJF_P99.9','P_day_MAM_P99.9')
  
#  P = calcPerformanceAttributes(clim=clim_ref,sim=sim,attSel=attSel,cSel = 1,vSel='P')
  #P = calcPerformanceAttributes(clim=clim_ref,sim=sim,attSel=attSel,cSel = 'mean',vSel='P')
P = calcPerformanceAttributes(clim=clim_ref,sim=sim,attSel=attSel)
  
  par(mfrow=c(5,5),mar=c(4,5,2,1))
  for (att in names(P)){
  #  plotPerformanceOAT(P, sim.addPET, metric=att,plotType='base',attSel=attPerturb)
    plotPerformanceOAT(P, sim, metric=att,plotType='base',attSel=attPerturb)
  } 
  
}

print(P)

pause

############################################################################

source(paste0(runDirname,'GR4J_funcs.R'))

dates = as.Date(clim_ref$times)
Param = setup_cal_GR4J(dates = dates,clim_ref_PET$P,clim_ref_PET$PET,Qobs)
systemArgs = list(dates=dates,Param=Param)
metrics = c('meanQ','P99','P25','min3yr')

sysOutSim = runSystemModel(sim=sim.addPET,systemModel = GR4J_wrapper,systemArgs = systemArgs,metrics = c('meanQ','P99','P25','min3yr'),varNames=c('P','PET'))
 
sysOutClim = GR4J_wrapper(data = clim_ref_PET, systemArgs = systemArgs,metrics=metrics)
  
eval = evaluate_system_metrics(sim=sim.addPET,clim=clim_ref_PET,
                        systemModel=GR4J_wrapper,systemArgs=systemArgs,
                        metrics=metrics,varNames=c('P','PET'))
par(mfrow=c(2,2))
for (metric in metrics){
  yAll = c(eval$systemPerf_base[[metric]],eval$systemPerf_obsClim[metric])
  ylim = c(min(yAll),max(yAll))
  boxplot.ext(as.vector(eval$systemPerf_base[[metric]]),ylim=ylim)
  points(eval$systemPerf_obsClim[metric],col='red')
  title(metric)
}

plotPerformanceOAT(performance = sysOutSim, sim=sim.addPET, metric = 'meanQ',attSel=attPerturb)
plotPerformanceOAT(performance = sysOutSim, sim=sim.addPET, metric = 'P99',attSel=attPerturb)
plotPerformanceOAT(performance = sysOutSim, sim=sim.addPET, metric = 'P25',attSel=attPerturb)
plotPerformanceOAT(performance = sysOutSim, sim=sim.addPET, metric = 'min3yr',attSel=attPerturb)

############################################################################

#dev.off()


pause

P=sim$Rep1$Target1$P$sim
clim_sim = clim_ref
#clim_sim$P = P[,1]
clim_sim$P = P[,1]

clim_ref1 = clim_ref
clim_ref1$P = clim_ref1$P[,1]


calculateAttributes(clim_sim,'P_day_all_P99.9') / calculateAttributes(clim_sim,'P_day_all_P99')
calculateAttributes(clim_ref1,'P_day_all_P99.9') / calculateAttributes(clim_ref1,'P_day_all_P99')
calculateAttributes(clim_sim,'P_day_SON_P99.9') / calculateAttributes(clim_sim,'P_day_SON_P99')
calculateAttributes(clim_ref1,'P_day_SON_P99.9') / calculateAttributes(clim_ref1,'P_day_SON_P99')
calculateAttributes(clim_sim,'P_day_DJF_P99.9') / calculateAttributes(clim_sim,'P_day_DJF_P99')
calculateAttributes(clim_ref1,'P_day_DJF_P99.9') / calculateAttributes(clim_ref1,'P_day_DJF_P99')

