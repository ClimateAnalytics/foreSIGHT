# script for generating single-site perturbed rainfall and PET data for Loxton for Source w crop model 

# clear all variables
rm(list=ls())

#----------------------------------------------------------------
# Set paths and get data
#----------------------------------------------------------------

# path to foreSIGHT code from repo  
foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT/'
#foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT_current_github/foreSIGHT/'

# load foreSIGHT package
setwd(foreSIGHTDir)
devtools::load_all()

# path to box folder
#boxDir = 'C:/Users/capta/Box/'
boxDir = 'C:/Users/a1065639/Box/'
# path for data files
dataDir <- paste0(boxDir,'2021 CRAFT Barossa project/')

# climate reference data file names
evapFile <- "External datasets/Hydroclimate (instrumental)/evap_1900to2020.Rdata"
rainFile <- "External datasets/Hydroclimate (instrumental)/rain_data_1900to2020_heggiesInfilled.Rdata"

# read in evap reference data
load(paste0(dataDir, evapFile))

evap.date = as.Date(evap_1900to2020[,1],format = '%d/%m/%Y')
evap.data = evap_1900to2020[,2]

# read in rain reference data
load(paste0(dataDir, rainFile))
nRain = dim(rain_data_1900to2020)[2]-1
rain.date = as.Date(rain_data_1900to2020[,1],format = '%d/%m/%Y')
rain.data = rain_data_1900to2020[,2:(nRain+1)]

if (!any(rain.date!=evap.date)){
  date=rain.date
} else {
  stop('dates for P and PET file differ')
}

# select time periods of data to use for simulations
startYr = 1976
endYr = 2005
# startYr = 1980
# endYr = 1989
catchment_keep = 1:3
year = as.integer(format(date,'%Y'))
keep = which((year>=startYr)&(year<=endYr))
# create reference data 
clim_ref <- list(times = as.POSIXct(date[keep],tz='UTC'),
                 P = rain.data[keep,catchment_keep],
                 PET = evap.data[keep])

######################################################################

modelSelection = list()

modelSelection$modelType = list()
modelSelection$modelType$P = "latent"

modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "seas"

modelSelection[["modelParameterBounds"]] <- list()
modelSelection[["modelParameterBounds"]][["P"]] <- list()
modelSelection[["modelParameterBounds"]][["P"]][["alpha.SON"]] <- c(0.1, 0.9)
modelSelection[["modelParameterBounds"]][["P"]][["alpha.DJF"]] <- c(0.1, 0.9)
modelSelection[["modelParameterBounds"]][["P"]][["alpha.MAM"]] <- c(0.1, 0.9)
modelSelection[["modelParameterBounds"]][["P"]][["alpha.JJA"]] <- c(0.1, 0.9)
modelSelection[["modelParameterBounds"]][["P"]][["sigma.SON"]] <- c(0.5, 8)
modelSelection[["modelParameterBounds"]][["P"]][["sigma.DJF"]] <- c(0.5, 8)
modelSelection[["modelParameterBounds"]][["P"]][["sigma.MAM"]] <- c(0.5, 8)
modelSelection[["modelParameterBounds"]][["P"]][["sigma.JJA"]] <- c(0.5, 8)
modelSelection[["modelParameterBounds"]][["P"]][["mu.SON"]] <- c(-8, 1)
modelSelection[["modelParameterBounds"]][["P"]][["mu.DJF"]] <- c(-8, 1)
modelSelection[["modelParameterBounds"]][["P"]][["mu.MAM"]] <- c(-8, 1)
modelSelection[["modelParameterBounds"]][["P"]][["mu.JJA"]] <- c(-8, 1)
modelSelection[["modelParameterBounds"]][["P"]][["lambda.SON"]] <- c(1, 5)
modelSelection[["modelParameterBounds"]][["P"]][["lambda.DJF"]] <- c(1, 5)
modelSelection[["modelParameterBounds"]][["P"]][["lambda.MAM"]] <- c(1, 5)
modelSelection[["modelParameterBounds"]][["P"]][["lambda.JJA"]] <- c(1, 5)


modelSelection[["optimisationArguments"]] <- list()
modelSelection[["optimisationArguments"]][["nMultiStart"]] <- 5
modelSelection[["optimisationArguments"]][["OFtol"]] <- 0.02
# modelSelection[["optimisationArguments"]][['RGN.control']]=list(iterMax=50)

modelSelection[["penaltyAttributes"]] <- c("P_day_all_tot", "P_day_all_xP99overPave",
                                          "P_day_all_nWet", "P_day_all_avgDSD")#,

modelSelection[["penaltyWeights"]] = rep(3,4)

modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")
write(modelSelectionJSON, file = controlFile)

#########################

attPerturb = c('P_day_all_avgDSD')
attHold = c('P_day_all_tot','P_day_all_xP99overPave','P_day_all_nWet')

attPerturbType = "regGrid"
attPerturbSamp = c(5)
attPerturbMin = c(1.)
attPerturbMax = c(1.4)

#########################

# attPerturb = c('P_day_all_tot')
# attHold = c('P_day_all_avgDSD','P_day_all_xP99overPave','P_day_all_nWet')
# 
# attPerturbType = "regGrid"
# attPerturbSamp = c(5)
# attPerturbMin = c(0.8)
# attPerturbMax = c(1.2)

#########################

# create the exposure space
expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType,
                          attHold = attHold)

attTied = list(seas='allTargets')

expSpace = tieAttributes(expSpace,attTied)

i=which(grepl(pattern = 'xP99overPave',x = colnames(expSpace$targetMat)))

expSpace$targetMat[,i] = 1/expSpace$targetMat[,'P_day_all_tot']

time.1 = Sys.time()
sim_stoch = generateScenarios(reference = clim_ref,
                              expSpace = expSpace,
                              controlFile = controlFile,
                              seedID = 1,
                              numReplicates = 1,
                              cores = 1)
time.2 = Sys.time()
print(time.2-time.1)

# fname = paste0('testing_2.0/Loxton_',attPerturb,'.RData') 
# save(file=fname,sim_stoch,clim_ref)

plotScenarios(sim_stoch)

attSel = colnames(expSpace$targetMat)

P = calcPerformanceAttributes(clim=clim_ref,sim=sim_stoch,attSel=attSel)

par(mfrow=c(5,5),mar=c(4,5,2,1))
for (att in names(P)){
  plotPerformanceOAT(P, sim_stoch, metric=att,plotType='base',attSel=attPerturb)
}

####################

clim_sim = list(times=clim_ref$times,P=sim_stoch$Rep5$Target1$P$sim)
calculateAttributes(clim_sim,'P_day_all_P99')
calculateAttributes(clim_ref,'P_day_all_P99')

####################

attSel = c('P_day_all_tot','P_day_all_avgDSD','P_day_all_P99','P_day_all_nWet')
P = calcPerformanceAttributes(clim=clim_ref,sim=sim_stoch,attSel=attSel)
for (att in names(P)){
  plotPerformanceOAT(P, sim_stoch, metric=att,plotType='base',attSel=attPerturb)
}
