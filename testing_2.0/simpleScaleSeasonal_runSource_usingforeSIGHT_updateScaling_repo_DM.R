# script for generating multi-site perturbed data using simple scaling with seasonality
# and use the multi-site data to run the eWater Source model
# perturbs original P and PET data and runs Source with this data

# clear all variables
rm(list=ls())

#----------------------------------------------------------------
# Set paths and get data
#----------------------------------------------------------------

# path to box folder
#boxDir = 'C:/Users/capta/Box/'
boxDir = 'C:/Users/a1065639/Box/'

# library(foreSIGHT)

# path to foreSIGHT code from repo  
#foreSIGHTDir = paste0(boxDir,'2021_foreSIGHT/git_current/MiddleRiverCRAFT_current/foreSIGHT/')
#foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT/'
foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT_current_github/foreSIGHT/'

# path for data files
dataDir <- paste0(boxDir,'2021 CRAFT Barossa project/')

# load foreSIGHT package
setwd(foreSIGHTDir)
devtools::load_all()

# climate reference data file names
evapFile <- "External datasets/Hydroclimate (instrumental)/evap_1900to2020.Rdata"
rainFile <- "External datasets/Hydroclimate (instrumental)/rain_data_1900to2020_heggiesInfilled.Rdata"

# read in evap reference data
load(paste0(dataDir, evapFile))

evap.date = as.Date(evap_1900to2020[,1],format = '%d/%m/%Y')
evap.data = evap_1900to2020[,2]

# evap_1900to2020 = as.vector(evap_1900to2020[,2])

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

#-----------------------------------------------------------------
# Step A: Create an exposure space
#-----------------------------------------------------------------

# specify perturbed attributes
attPerturb <- c("PET_day_all_avg_m", "P_day_all_tot_m","P_day_all_seasRatioMarMay") # note that updates to attribute manager allow customisation of months used for seasRatio 

# specify perturbation type and minimum-maximum ranges of the perturbed attributes
attPerturbType <- "regGrid"
#attPerturbSamp <- c(3,5,5)     
#attPerturbMin <- c(1, 0.7,0.7)
#attPerturbMax <- c(1.2, 1.3,1.3)
attPerturbSamp <- c(2,2,2)     
attPerturbMin <- c(1, 0.7,0.7)
attPerturbMax <- c(1.2, 1.3,1.3)
# create the exposure space using foreSIGHT
expSpace <- createExpSpace(attPerturb = attPerturb,
                           attPerturbSamp = attPerturbSamp,            # set to null
                           attPerturbMin = attPerturbMin,
                           attPerturbMax = attPerturbMax,
                           attPerturbType = attPerturbType,
                           attHold = NULL)


#-----------------------------------------------------------------
# Step B: Generate perturbed time series
#-----------------------------------------------------------------

# select time periods of data to use for simulations
startYr = 1976
endYr = 2005
year = as.integer(format(date,'%Y'))
keep = which((year>=startYr)&(year<=endYr))
# create reference data 
clim_ref <- list(times = as.POSIXct(date[keep],tz='UTC'),
                 P = rain.data[keep,],
                 PET = evap.data[keep])

# generate simulations using simple scaling with seasonality 
sim <- generateScenarios(reference = clim_ref,             # input observed data
                         expSpace = expSpace,        # exposure space created by the user
                         controlFile = "scaling")    # using simple scaling with seasonality

#save(sim,file=paste0(dataDir,"SWModel/BarossaModelFiles/Outputs/1.4_SeasonalityChange_Run/sim.Rdata"))

# # do some checks
# month = as.integer(format(clim_ref$times,'%m'))
# keep.seas = which(month>=3 & month<=5)
# 
# attObs = calculateAttributes(climateData = clim_ref,attSel = attPerturb)
# 
# P = clim_ref$P
# Pseas = apply(P[keep.seas,],2,sum)
# Pall = apply(P,2,sum)
# seas_ratio_obs = Pseas/Pall
# 
# P = sim$Rep1$Target1$P
# Pseas = apply(P[keep.seas,],2,sum)
# Pall = apply(P,2,sum)
# seas_ratio_sim = Pseas/Pall
# 
# seas_ratio_sim/seas_ratio_obs
#   
# #date = as.Date(paste0(clim_ref$year,'/',clim_ref$month,'/',clim_ref$day))
# date = clim_ref$times
# dateOneYear = date[1:365]
# doy = as.integer(format(date,'%j'))
# clim = calc_ClimDaily_dayOfYearWindow(obs=clim_ref$P[,1],dateObs = date,dateClim=dateOneYear,inc = 14)
# meanClim = apply(clim,1,mean)
# plot(meanClim,col='red',type='l',ylim=c(0,3.5))
# 
# for (target in names(sim$Rep1)){
#   print(target)
#   P=sim$Rep1[[target]]$P
#   PET=sim$Rep1[[target]]$PET
#   simData = list(times=clim_ref$times,
#                  P=P,
#                  PET=PET)
#   
#   attSim = calculateAttributes(climateData = simData,attSel = attPerturb)
#   print(attSim/attObs)
#   
#   Pseas = apply(P[keep.seas,],2,sum)
#   Pall = apply(P[,],2,sum)
#   seas_ratio_sim = Pseas/Pall
#   print(seas_ratio_sim/seas_ratio_obs)
#   
#   clim = calc_ClimDaily_dayOfYearWindow(obs=P[,1],dateObs = date,dateClim=dateOneYear,inc = 14)
#   meanClim = apply(clim,1,mean)
#   lines(meanClim,col='gray',lwd=0.6)
#   
# }

#-----------------------------------------------------------------
# Step C: Run the system model
#-----------------------------------------------------------------

#DM: I have not changed code from here onwards

# load the source wrapper function
source(paste0(boxDir, "2021 CRAFT Barossa project/Code/source_wrapper.R"))

# directory containing source input datasets
inputDir <- paste0(boxDir, "2021 CRAFT Barossa project/SWModel/BarossaModelFiles/Inputs/")

# names of climate data files to be overwritten; perturbed time series will be written to these files
# these are the files that Source reloads on run, so their names have to be the same
evapOutFile <- "Input_evap_Original.csv"
rainOutFile <- "Input_rain_Original.csv"

# source variables
source_exe_path <- "C:/Program Files/eWater/Source 5.0.0.10962/"
source_exe <- "RiverSystem.CommandLine.exe"
project_path <- paste0(boxDir, "2021 CRAFT Barossa project/SWModel/BarossaModelFiles/Scenarios/Barossa_Sce2_stressTest1_V5.0.rsproj")
scenario <- "Scenario 1"
runconfiguration <-"Single analysis"
output_file <- paste0(boxDir, "2021 CRAFT Barossa project/SWModel/BarossaModelFiles/Outputs/BaselineRun_Aug17/source_output") #changed in runSystemModelFunction
runstarttime <- "1/1/1976"
runendtime <- "31/12/2005"

#can create new functions for other types such as max - specify in systemArgs$outSummaryFunc below
meanAnnual <- function(x, date) {
  annAgg <- aggregate(x, by = list(as.POSIXlt(date)$year), FUN=sum)[["x"]]  # annual totals
  annAggMean <- mean(annAgg)
  return(annAggMean)
}

systemArgs <- list()
systemArgs$PFile <- paste0(inputDir, rainOutFile)
systemArgs$PETFile <- paste0(inputDir, evapOutFile)
systemArgs$startDate <- paste0(startYr, "/01/01")
systemArgs$endDate <- paste0(endYr, "/12/31")
systemArgs$P_colNames <- colnames(rain_data_1900to2020)
systemArgs$PET_colNames <- c("Date", "PET")
systemArgs$source_exe_path <- source_exe_path
systemArgs$source_exe <- source_exe
systemArgs$project_path <- project_path
systemArgs$scenario <- scenario
systemArgs$runconfiguration <- runconfiguration
systemArgs$runstarttime <- runstarttime
systemArgs$runendtime <- runendtime
systemArgs$outFilename <- output_file
# specify chosen Source results based on column named when imported to R
systemArgs$outVarname <- c("Downstream.Flow", "Supplied.Demand.Volume", "Actual.Demand.Volume", "Actual.Evaporation.Volume", "Storage.Record", "Total.Water.Supplied", "Volume.Ordered")
systemArgs$outSummaryFunc <- meanAnnual

metricNames <- c("Annual Yaldara Flow (ML)", "Annual Farm Dam Demand Supplied (ML)", "Annual Farm Dam Demand (ML)", "Annual Farm Dam Evaporation (ML)", "Annual Farm Dam Storage (ML)",
                 "Annual Water Course Demand Supplied (ML)", "Annual Water Course Demand (ML)")

#################################

# To test source_wrapper if needed
testOutput <- source_wrapper(data = sim$Rep1$Target1,
                             systemArgs = systemArgs)

#################################

# allow output filename to be change to include rep and target info
systemArgs$outFilename_tarRepStr='long'

systemPerf = runSystemModel(sim = sim,              # simulation; the perturbed time series
                            systemModel = source_wrapper,      # the system model function
                            systemArgs = systemArgs,           # argument to the system model function
                            metrics = metricNames)             # selected performance metrics 


#################################

# runSystemModelMultipleOut <- function (sim, systemModel, systemArgs, metrics){
#   origOutName <- systemArgs$outFilename
#   repNames <- names(sim[grep("Rep", names(sim))])
#   tarNames <- names(sim[[repNames[1]]])
#   nRep <- length(repNames)
#   nTar <- length(tarNames)
#   temp <- names(sim[[repNames[1]]][[tarNames[1]]])
#   varNames <- temp #temp[which(temp %in% fSVars)]
#   rm(temp)
#   performance <- vector("list", length = length(metrics))
#   for (i in 1:length(metrics)) performance[[i]] <- matrix(NA, 
#                                                           nrow = nTar, ncol = nRep)
#   for (r in 1:nRep) {
#     for (t in 1:nTar) {
#       scenarioData <- sim[["simDates"]]
#       varTemp <- list()
#       for (v in varNames) {
#         if ((is.character(sim[["controlFile"]]))) {
#           if (sim[["controlFile"]] == "scaling") {
#             max_nSites <- max(sapply(sim[[repNames[r]]][[tarNames[t]]], 
#                                      ncol))
#             if (max_nSites == 1) {
#               varTemp <- as.data.frame(sim[[repNames[r]]][[tarNames[t]]][[v]])
#             }
#             else {
#               varTemp[[v]] <- sim[[repNames[r]]][[tarNames[t]]][[v]]
#             }
#           }
#           else {
#             stop(paste0("sim$controlFile unrecognized."))
#           }
#         }
#         else {
#           max_nSites <- max(sapply(sim[[repNames[r]]][[tarNames[t]]], 
#                                    function(x) {
#                                      ncol(x[["sim"]])
#                                    }))
#           if (max_nSites == 1) {
#             varTemp <- as.data.frame(sim[[repNames[r]]][[tarNames[t]]][[v]][["sim"]])
#           }
#           else {
#             varTemp[[v]] <- sim[[repNames[r]]][[tarNames[t]]][[v]][["sim"]]
#           }
#         }
#         if (is.data.frame(varTemp)) {
#           names(varTemp) <- v
#           scenarioData <- cbind(scenarioData, varTemp)
#           rm(varTemp)
#         }
#         else {
#           scenarioData <- varTemp
#         }
#       }
#       outname <- paste0(origOutName,repNames[r],"_",tarNames[t],"_",varNames[1],sim$expSpace$targetMat[t,1],"_",varNames[2],sim$expSpace$targetMat[t,2],".csv")
#       systemArgs$outFilename <- outname
#       perfTemp <- systemModel(data = scenarioData, systemArgs = systemArgs, 
#                               metrics = metrics)
#       for (i in 1:length(metrics)) performance[[i]][t, 
#                                                     r] <- perfTemp[[i]]
#       rm(scenarioData)
#     }
#   }
#   names(performance) <- metrics
#   return(performance)
# }
# 
# systemPerf <- runSystemModelMultipleOut(sim = sim,              # simulation; the perturbed time series
#                                         systemModel = source_wrapper,      # the system model function
#                                         systemArgs = systemArgs,           # argument to the system model function
#                                         metrics = metricNames)             # selected performance metrics 

#--------------------------------------------------------------
# Step D: Plot the performance metrics
#--------------------------------------------------------------

simSummary <- getSimSummary(sim)

p1 <- plotPerformanceSpace(systemPerf[1], simSummary)
p2 <- plotPerformanceSpace(systemPerf[2], simSummary)
p3 <- plotPerformanceSpace(systemPerf[3], simSummary)
p4 <- plotPerformanceSpace(systemPerf[4], simSummary)

# Coverting the sum of farm dam storage to average
systemPerfTemp=systemPerf[5]
systemPerfTemp[[1]]=systemPerf[5][[1]]/365.25
p5 <- plotPerformanceSpace(systemPerfTemp, simSummary)

p6 <- plotPerformanceSpace(systemPerf[6], simSummary)
p7 <- plotPerformanceSpace(systemPerf[7], simSummary)


p1_oat <- plotPerformanceOAT(systemPerf[1], simSummary)
p2_oat <- plotPerformanceOAT(systemPerf[2], simSummary)
p3_oat <- plotPerformanceOAT(systemPerf[3], simSummary)
p4_oat <- plotPerformanceOAT(systemPerf[4], simSummary)
p5_oat <- plotPerformanceOAT(systemPerfTemp, simSummary)
p6_oat <- plotPerformanceOAT(systemPerf[6], simSummary)
p7_oat <- plotPerformanceOAT(systemPerf[7], simSummary)



