rm(list=ls())

#foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT/'
foreSIGHTDir = '/scratchdata1/users/a1065639/DEW_foreSIGHT/foreSIGHT/'

#devtools::load_all(foreSIGHTDir)
library(foreSIGHT)

runDirname = paste0(foreSIGHTDir,'testing_2.0/')
setwd(runDirname)

source(paste0(runDirname,'GR4J_funcs.R'))
source(paste0(runDirname,'boxplot.ext_DM.r'))

############################################################################

catchment = 'A5030502' # Scott Creek
#catchment = 'A5050517' # North Para River at Penrice (A5050517)
startYr = 1976
endYr = 1985
#endYr = 1995
#endYr = 2005

numReplicates = 50
cores = 50

load_data = TRUE # read RData files (TRUE) or create them using DroughtRisk Package (FALSE) 

############################################################################

runStr = 'V1'
#runStr = paste0(runStr,'_',catchment,'_',startYr,'_',endYr,
#		paste(attPerturb,collapse='_'))

#fname = paste0(runDirname,'summary_',catchment,'_',startYr,'_',endYr,'.pdf')
#fname = paste0(runDirname,'sim_',catchment,'_',startYr,'_',endYr,
#	       paste(attPerturb,collapse='_'),'.RData')


############################################################################

fname = paste0(runDirname,'data_',catchment,'_',startYr,'_',endYr,'.RData')
if (!load_data){
  source('load_camels.R')
  devtools::load_all('/Users/a1065639/Work/DroughtRisk/')
  data = load_camels(catchment,startYr,endYr,climPET = T)
  save(file=fname,data)
} else {
  load(file=fname)
}

# create reference data 
clim_ref <- list(times = data$times,
                 P = data$P)  

############################################################################

#fname = paste0(runDirname,'summary_',catchment,'_',startYr,'_',endYr,'.pdf')
#fname = paste0(runDirname,'summary_',runStr,'.pdf')
#pdf(fname)

############################################################################

attPerturbType = "regGrid"
attPerturb = c('P_day_all_seasRatioMarMay','P_day_all_P99')
#attPerturbSamp = c(1,1)
#attPerturbMin = c(1,1)
#attPerturbMax = c(1,1)
attPerturbSamp = c(5,5)
attPerturbMin = c(0.7,1.)
attPerturbMax = c(1.3,1.3)
attHold = c('P_day_all_tot','P_day_all_avgDSD','P_day_all_nWet',
            'P_day_DJF_avgDSD','P_day_MAM_avgDSD','P_day_JJA_avgDSD','P_day_SON_avgDSD',
            'P_day_DJF_nWet','P_day_MAM_nWet','P_day_JJA_nWet','P_day_SON_nWet',
            'P_day_DJF_tot','P_day_JJA_tot','P_day_SON_tot')
attTied = list(P_day_all_P99=c('P_day_DJF_normP99','P_day_MAM_normP99','P_day_JJA_normP99','P_day_SON_normP99'),
              P_day_all_seasRatioMarMay=c('P_day_MAM_tot'))

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
modelSelection$modelType$P = "latent"

modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "seas"

modelSelection[["optimisationArguments"]] <- list()
modelSelection[["optimisationArguments"]][["OFtol"]] <- 0.1

modelSelection[["penaltyAttributes"]] <- c('P_day_all_seasRatioMarMay','P_day_all_P99','P_day_all_tot',
                                           'P_day_all_avgDSD','P_day_all_nWet')
modelSelection[["penaltyWeights"]] = c(3,3,3,1.5,1.5)

modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")
write(modelSelectionJSON, file = controlFile)

############################################################################

runStr = paste0(runStr,'_',catchment,'_',startYr,'_',endYr,
                paste(attPerturb,collapse='_'))

#fname = paste0(runDirname,'summary_',catchment,'_',startYr,'_',endYr,'.pdf')
fname = paste0(runDirname,'summary_',runStr,'.pdf')
pdf(fname)

############################################################################

#time.1 = Sys.time()
#sim = generateScenarios(reference = clim_ref,
#                        expSpace = expSpace,
#                        controlFile = controlFile,
#                        seedID = 1,
#                        numReplicates = numReplicates,
#                        cores = cores)
#time.2 = Sys.time()
#print(time.2-time.1)

#fname = paste0(runDirname,'summary_',paste(attPerturb,collapse='_'),'.RData')
fname = paste0(runDirname,'sim_',runStr,'.RData')


#save.image(file=fname)

load(fname)

##########################################################################

plotScenarios(sim)

############################################################################

nTar = length(sim$Rep1)
if (nTar>1){
  
  attSel = colnames(expSpace$targetMat)
  attSel = c(attSel,'P_day_all_tot_dwellTime','P_day_all_avgWSD',
             'P_day_all_P99.9','P_day_JJA_P99.9','P_day_SON_P99.9',
             'P_day_DJF_P99.9','P_day_MAM_P99.9')
  
  for (att in attPerturb){
    par(mfrow=c(5,5),mar=c(4,5,2,1))
    plotPerformanceAttributesOAT(clim=clim_ref,
                                 sim=sim,
                                 attPerturb=att,
                                 attEval=attSel)    
  }

}

############################################################################


Qobs = data$Qobs
PET = data$PET 

dates = as.Date(clim_ref$times)
Param = setup_cal_GR4J(dates = dates,P=clim_ref$P,PET=PET,Qobs=Qobs)
systemArgs = list(dates=dates,Param=Param,PET=PET)
metrics = c('meanQ','P99','P25','min3yr')

sysOutSim = runSystemModel(sim=sim,systemModel = GR4J_wrapper,systemArgs = systemArgs,metrics = metrics)
 
# print('calc sysOutClim')
# sysOutClim = GR4J_wrapper(data = clim_ref, systemArgs = systemArgs,metrics=metrics)
# print('done')
# 
print('calc evaluate_system_metrics')
eval = evaluate_system_metrics(sim=sim,clim=clim_ref,
                               systemModel=GR4J_wrapper,systemArgs=systemArgs,
                               metrics=metrics)

par(mfrow=c(2,2))
for (metric in metrics){
  yAll = c(eval$systemPerf_base[[metric]],eval$systemPerf_obsClim[metric])
  ylim = c(min(yAll),max(yAll))
  boxplot.ext(as.vector(eval$systemPerf_base[[metric]]),ylim=ylim)
  points(eval$systemPerf_obsClim[metric],col='red')
  title(metric)
}

if (nTar>1){
  
  plotPerformanceOAT(performance = sysOutSim, sim=sim, metric = 'meanQ',attSel=attPerturb[1])
  plotPerformanceOAT(performance = sysOutSim, sim=sim, metric = 'P99',attSel=attPerturb[1])
  plotPerformanceOAT(performance = sysOutSim, sim=sim, metric = 'P25',attSel=attPerturb[1])
  plotPerformanceOAT(performance = sysOutSim, sim=sim, metric = 'min3yr',attSel=attPerturb[1])

  plotPerformanceOAT(performance = sysOutSim, sim=sim, metric = 'meanQ',attSel=attPerturb[2])
  plotPerformanceOAT(performance = sysOutSim, sim=sim, metric = 'P99',attSel=attPerturb[2])
  plotPerformanceOAT(performance = sysOutSim, sim=sim, metric = 'P25',attSel=attPerturb[2])
  plotPerformanceOAT(performance = sysOutSim, sim=sim, metric = 'min3yr',attSel=attPerturb[2])

  plotPerformanceSpace(performance = sysOutSim, sim=sim, metric = 'meanQ',type='filled.contour',nContour=5)
  plotPerformanceSpace(performance = sysOutSim, sim=sim, metric = 'P99',type='filled.contour',nContour=5)
  plotPerformanceSpace(performance = sysOutSim, sim=sim, metric = 'P25',type='filled.contour',nContour=5)
  plotPerformanceSpace(performance = sysOutSim, sim=sim, metric = 'min3yr',type='filled.contour',nContour=5)
  
}

############################################################################

dev.off()



