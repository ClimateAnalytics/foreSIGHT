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
#library(foreSIGHT)

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
year = as.integer(format(date,'%Y'))
keep = which((year>=startYr)&(year<=endYr))
# create reference data 
clim_ref <- list(times = as.POSIXct(date[keep],tz='UTC'),
                 P = P[keep],
                 PET = PET[keep])

###################################

attPerturb = c('P_day_all_normP99',"PET_day_all_avg")
attHold = c('P_day_all_tot_m','P_day_all_avgDSD','P_day_all_nWet')

attPerturbType = "regGrid"
#attPerturbSamp = c(2,2)
#attPerturbMin = c(1,1)
#attPerturbMax = c(1.3,1.3)
attPerturbSamp = c(1,1)
attPerturbMin = c(1.3,1.3)
attPerturbMax = c(1.3,1.3)

expSpace.1 = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType,
                          attHold = attHold)


attsTied = setSeasonalTiedAttributes(c('P_day_all_normP99','P_day_all_tot_m',
                                       'P_day_all_avgDSD','P_day_all_nWet'))

expSpace.2 = tieAttributes(expSpace.1,attsTied)

attsTied = list(PET_day_all_avg=c("mv.PET.P_day_all_avgDryDay",
                                  "mv.PET.P_day_all_avgWetDay"))

expSpace = tieAttributes(expSpace.2,attsTied)

###################################

modelSelection = list()

modelSelection$modelType = list()
modelSelection$modelType$P = "latent"
modelSelection$modelType$PET = "wgenO"

modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "seas"
modelSelection$modelParameterVariation$PET = "harWD"

o = modCalibrator(obs=clim_ref,modelTag='PET-harWD-wgenO')

modelSelection[["modelParameterBounds"]] <- list()
modelSelection[["modelParameterBounds"]][["PET"]] <- list()

parPET = o[[1]]
parPET = parPET[!names(parPET)%in%c('mu.D.m','mu.W.m')]
for (par in names(parPET)){
  modelSelection[["modelParameterBounds"]][["PET"]][[par]] <- c(parPET[par], parPET[par])
}

modelSelection[["optimisationArguments"]] <- list()
modelSelection[["optimisationArguments"]][["nMultiStart"]] <- 5
modelSelection[["optimisationArguments"]][["OFtol"]] <- 0.02

modelSelection[["penaltyAttributes"]] <- c("P_day_all_normP99","P_day_all_tot_m", 
                                           "P_day_all_nWet", "P_day_all_avgDSD",
                                           "PET_day_all_avg")
modelSelection[["penaltyWeights"]] <-c(5,5,2,2,2)

modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")
write(modelSelectionJSON, file = controlFile)

####################

time.1 = Sys.time()
sim_stoch = generateScenarios(reference = clim_ref,
                              expSpace = expSpace,
                              controlFile = controlFile,
                              seedID = 1,
                              numReplicates = 3,
                              cores = 1)
time.2 = Sys.time()
print(time.2-time.1)

fname = paste0('testing_2.0/Loxton_V2_',attPerturb,'.RData') 
save(file=fname,sim_stoch,clim_ref)

plotScenarios(sim_stoch)

attSel = colnames(expSpace$targetMat)

pause





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


####################

