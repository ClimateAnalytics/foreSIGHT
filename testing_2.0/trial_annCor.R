rm(list=ls())

#####################################################################
# setup data

devtools::load_all('/Users/a1065639/Work/DroughtRisk/')

#catchment = '410057'
catchment = 'A5030502'

# startYr = 1976
# endYr = 2005

# startYr = 1920
# endYr = 2010

startYr = 1976
endYr = 1985

data = DroughtRisk::read_catchment_data(catchment=catchment,dataSource = 'CAMELS-AUS')

year = as.integer(format(data$dates,'%Y'))
keep = which((year>=startYr)&(year<=endYr))

dates = data$dates[keep]
P = data$P[keep]
PET = data$PETobs[keep]
Qobs = data$Qobs[keep]

# create reference data 
clim_ref <- list(times = as.POSIXct(dates,tz='UTC'),
                 P = P,
                 PET = PET)

#####################################################################

# path to foreSIGHT code from repo  
foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT/'
#foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT_current_github/foreSIGHT/'

# path for data files
#dataDir = 'C:/Users/a1065639/Box/2025_DEW_foreSIGHT/Data/'

# load foreSIGHT package
#setwd(foreSIGHTDir)
devtools::load_all()
#devtools::install()
#library(foreSIGHT)



######################################################################

modelSelection = list()

modelSelection$modelType = list()
modelSelection$modelType$P = "latent"

modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "seas"

modelSelection[["optimisationArguments"]] <- list()
modelSelection[["optimisationArguments"]][["nMultiStart"]] <- 5
modelSelection[["optimisationArguments"]][['RGN.control']]=list(iterMax=50)

modelSelection[["penaltyAttributes"]] <- c("P_day_all_tot", "P_day_all_xP99overPave",
                                           "P_day_all_nWet", "P_day_all_avgDSD")#,"P_year_all_cv")
#modelSelection[["penaltyWeights"]] = rep(2,5)
modelSelection[["penaltyWeights"]] = rep(2,4)

modelSelection$postProcessing = list(P=list())
#modelSelection$postProcessing$P$types = c('annVar')
#modelSelection$postProcessing$P$types = c('scaleExtremesSeas')

modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")
write(modelSelectionJSON, file = controlFile)

####################

attPerturb = c('P_day_all_xP99overPave')
attHold = c('P_day_all_tot','P_day_all_avgDSD','P_day_all_nWet')#,'P_year_all_cv')

attPerturbType = "regGrid"
attPerturbSamp = c(5)
attPerturbMin = c(0.9)
attPerturbMax = c(1.3)

#########################

targetTypes = NULL

# create the exposure space
expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType,
                          attHold = attHold,
                          targetTypes = targetTypes)

#attTied = list(seas='allTargets')
attTied = list(seas=c('P_day_all_tot','P_day_all_avgDSD',
                      'P_day_all_xP99overPave','P_day_all_nWet'))

expSpace = tieAttributes(expSpace,attTied)

# i=which(grepl(pattern = 'xP99overPave',x = colnames(expSpace$targetMat)))
# expSpace$targetMat[,i] = 1/expSpace$targetMat[,'P_day_all_tot']


time.1 = Sys.time()
sim_stoch = generateScenarios(reference = clim_ref,
                              expSpace = expSpace,
                              controlFile = controlFile,
                              seedID = 1,
                              numReplicates = 1)
time.2 = Sys.time()
print(time.2-time.1)

plotScenarios(sim_stoch)


attSel = colnames(sim_stoch$expSpace$targetMat)
par(mfrow=c(4,4))
P = calcPerformanceAttributes(clim=clim_ref,sim=sim_stoch,attSel=attSel)
for (att in names(P)){
  plotPerformanceOAT(P, sim_stoch, metric=att,plotType='base',attSel=attPerturb)
}


#########################

