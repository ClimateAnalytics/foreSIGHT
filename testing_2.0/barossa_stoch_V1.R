
rm(list=ls())

foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT/'
#foreSIGHTDir = '/scratchdata1/users/a1065639/DEW_foreSIGHT/foreSIGHT/'

devtools::load_all(foreSIGHTDir)
#devtools::install(foreSIGHTDir)

runDirname = paste0(foreSIGHTDir,'testing_2.0/')
setwd(runDirname)

############################################################################

#load('data_A5050517_1976_2005.RData')
#numReplicates = 20
#cores = 25

numReplicates = 1
cores = 1

#clim_ref = convert_climYMD_POSIXct(barossa_obs)
#clim_ref = barossa_obs
#clim_ref$P = clim_ref$P[,1:3]

load('barossa_clim_ref.RData')

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
modelSelection[["optimisationArguments"]][["nMultiStart"]] <- 1
modelSelection[["optimisationArguments"]][["OFtol"]] <- 0.05
#modelSelection[["optimisationArguments"]][['RGN.control']] = list(iterMax=100)
modelSelection[["optimisationArguments"]][['RGN.control']] = list(iterMax=1)

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
#attHold = c('P_day_all_tot','P_day_all_P99','P_day_all_nWet','P_day_all_tot_cv')

attPerturbType = "regGrid"
attPerturb = c('P_day_all_tot')
attPerturbSamp = c(3)
attPerturbMin = c(0.8)
attPerturbMax = c(1.2)
attHold = c('P_day_all_avgDSD','P_day_all_P99','P_day_all_nWet','P_day_all_tot_cv')

#attPerturbType = "regGrid"
#attPerturb = c('P_day_all_P99')
#attPerturbSamp = c(5)
#attPerturbMin = c(0.9)
#attPerturbMax = c(1.3)
#attHold = c('P_day_all_tot','P_day_all_avgDSD','P_day_all_nWet','P_day_all_tot_cv')

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

pdf(file=paste0(runDirname,'barossa_summary_',attPerturb,'.pdf'))

#plotScenarios(sim)

############################################################################

nTar = length(sim$Rep1)
if (nTar>1){

  attSel = colnames(expSpace$targetMat)
  attSel = c(attSel,'P_day_all_tot_dwellTime','P_day_all_P99.9','P_day_all_avgWSD')

  P = calcPerformanceAttributes(clim=clim_ref,sim=sim,attSel=attSel,cSel='mean',vSel='P')

  par(mfrow=c(5,5),mar=c(4,5,2,1))
  for (att in names(P)){
    plotPerformanceOAT(P, sim, metric=att,plotType='base',attSel=attPerturb)
  }

}

##########################################################################

dev.off()

save.image(file=paste0(runDirname,'barossa_summary_',attPerturb,'.RData'))


