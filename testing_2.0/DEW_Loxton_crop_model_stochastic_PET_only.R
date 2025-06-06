# script for generating single-site perturbed rainfall and PET data for Loxton for Source w crop model 

# clear all variables
rm(list=ls())

#----------------------------------------------------------------
# Set paths and get data
#----------------------------------------------------------------

# path to foreSIGHT code from repo  
foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT/'
#foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT_current_github/foreSIGHT/'

# path for data files
dataDir = 'C:/Users/a1065639/Box/2025_DEW_foreSIGHT/Data/'

# load foreSIGHT package
setwd(foreSIGHTDir)
devtools::load_all()

# remove.packages("foreSIGHT")
# .rs.restartR()
# devtools::install()
library(foreSIGHT)

# climate reference data file name
fname = paste0(dataDir,'24007_LOXTON_Dec2023.csv')

# read in evap reference data
d = read.csv(fname,header=T,skip=4017) # note there is different format for dates pre and post 1900. have chosen to not read in pre-1900 data.

date = as.Date(d[,1],format = '%d/%m/%Y')
P = d[,2]
PET = d[,3]

# select time periods of data to use for simulations
startYr = 1976
endYr = 2005
# startYr = 1980
# endYr = 1989
year = as.integer(format(date,'%Y'))
keep = which((year>=startYr)&(year<=endYr))
# create reference data 
clim_ref <- list(times = as.POSIXct(date[keep],tz='UTC'),
                 P = P[keep],
                 PET = PET[keep])

#-----------------------------------------------------------------

######################################################################


modelSelection = list()

modelSelection$modelType = list()
modelSelection$modelType$PET = "wgenLM"

modelSelection$modelParameterVariation = list()
#modelSelection$modelParameterVariation$PET = "seas"
modelSelection$modelParameterVariation$PET = "har"


# modelSelection[["modelParameterBounds"]][["PET"]] <- list()
# modelSelection[["modelParameterBounds"]][["PET"]][["sigma.SON"]] <- c(0.01, 5)
# modelSelection[["modelParameterBounds"]][["PET"]][["sigma.DJF"]] <- c(0.01, 5)
# modelSelection[["modelParameterBounds"]][["PET"]][["sigma.MAM"]] <- c(0.01, 5)
# modelSelection[["modelParameterBounds"]][["PET"]][["sigma.JJA"]] <- c(0.01, 5)

# modelSelection[["modelParameterBounds"]][["PET"]][["cor0.SON"]] <- c(0.1, 0.99)
# modelSelection[["modelParameterBounds"]][["PET"]][["cor0.DJF"]] <- c(0.1, 0.99)
# modelSelection[["modelParameterBounds"]][["PET"]][["cor0.MAM"]] <- c(0.1, 0.99)
# modelSelection[["modelParameterBounds"]][["PET"]][["cor0.JJA"]] <- c(0.1, 0.99)
# modelSelection[["modelParameterBounds"]][["PET"]][["muW.SON"]] <- c(0, 10)
# modelSelection[["modelParameterBounds"]][["PET"]][["muW.DJF"]] <- c(0, 10)
# modelSelection[["modelParameterBounds"]][["PET"]][["muW.MAM"]] <- c(0, 10)
# modelSelection[["modelParameterBounds"]][["PET"]][["muW.JJA"]] <- c(0, 10)
# modelSelection[["modelParameterBounds"]][["PET"]][["muD.SON"]] <- c(0, 15)
# modelSelection[["modelParameterBounds"]][["PET"]][["muD.DJF"]] <- c(0, 15)
# modelSelection[["modelParameterBounds"]][["PET"]][["muD.MAM"]] <- c(0, 15)
# modelSelection[["modelParameterBounds"]][["PET"]][["muD.JJA"]] <- c(0, 15)
# modelSelection[["modelParameterBounds"]][["PET"]][["sigmaD.SON"]] <- c(0.01, 5)
# modelSelection[["modelParameterBounds"]][["PET"]][["sigmaD.DJF"]] <- c(0.01, 5)
# modelSelection[["modelParameterBounds"]][["PET"]][["sigmaD.MAM"]] <- c(0.01, 5)
# modelSelection[["modelParameterBounds"]][["PET"]][["sigmaD.JJA"]] <- c(0.01, 5)


modelSelection[["optimisationArguments"]] <- list()
modelSelection[["optimisationArguments"]][["nMultiStart"]] <- 5
# modelSelection[["optimisationArguments"]][['RGN.control']]=list(iterMax=50)

modelSelection[["penaltyAttributes"]] <- c("PET_day_all_avg","PET_day_all_cor","PET_day_all_cv")
modelSelection[["penaltyWeights"]] = rep(3,3)
# modelSelection[["penaltyWeights"]] = rep(0,3)


modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")
write(modelSelectionJSON, file = controlFile)

####################

func_cv = function(data){
  return(sd(data,na.rm=T)/mean(data,na.rm=T))
}

#########################

attPerturb = c("PET_day_all_avg")
attHold = c("PET_day_all_cor","PET_day_all_cv")

attPerturbType = "regGrid"
attPerturbSamp = c(5)
attPerturbMin = c(1)
attPerturbMax = c(1.4)

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




time.1 = Sys.time()
sim_stoch = generateScenarios(reference = clim_ref,
                              expSpace = expSpace,
                              controlFile = controlFile,
                              seedID = 1,
                              numReplicates = 5,
                              cores = 1)
time.2 = Sys.time()
print(time.2-time.1)

plotScenarios(sim_stoch)

attSel = colnames(expSpace$targetMat)

P = calcPerformanceAttributes(clim=clim_ref,sim=sim_stoch,attSel=attSel)

par(mfrow=c(5,5),mar=c(4,5,2,1))
for (att in names(P)){
  plotPerformanceOAT(P, sim_stoch, metric=att,plotType='base',attSel=attPerturb)
}

####################

