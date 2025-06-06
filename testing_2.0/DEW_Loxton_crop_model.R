# script for generating single-site perturbed rainfall and PET data for Loxton for Source w crop model 

# clear all variables
rm(list=ls())

#----------------------------------------------------------------
# Set paths and get data
#----------------------------------------------------------------

# path to foreSIGHT code from repo  
#foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT/'
foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT_current_github/foreSIGHT/'

# path for data files
dataDir = 'C:/Users/a1065639/Box/2025_DEW_foreSIGHT/Data/'

# load foreSIGHT package
setwd(foreSIGHTDir)
devtools::load_all()

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

#-----------------------------------------------------------------
# Step A: Create an exposure space
#-----------------------------------------------------------------

# specify perturbed attributes
attPerturb <- c("PET_day_all_avg_m", "P_day_all_tot_m") 

# specify perturbation type and minimum-maximum ranges of the perturbed attributes
attPerturbType <- "regGrid"
attPerturbSamp <- c(7,7)     
attPerturbMin <- c(1, 0.7)
attPerturbMax <- c(1.3, 1.3)
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



# generate simulations using simple scaling with seasonality 
sim <- generateScenarios(reference = clim_ref,             # input observed data
                         expSpace = expSpace,        # exposure space created by the user
                         controlFile = "scaling")    # using simple scaling with seasonality



######################################################################



######################################################################

# #-----------------------------------------------------------------
# # Step C: Run the system model
# #-----------------------------------------------------------------
# 
# #DM: I have not changed code from here onwards
# 
# # load the source wrapper function
# source(paste0(boxDir, "2021 CRAFT Barossa project/Code/source_wrapper.R"))
# 
# # directory containing source input datasets
# inputDir <- paste0(boxDir, "2021 CRAFT Barossa project/SWModel/BarossaModelFiles/Inputs/")
# 
# # names of climate data files to be overwritten; perturbed time series will be written to these files
# # these are the files that Source reloads on run, so their names have to be the same
# evapOutFile <- "Input_evap_Original.csv"
# rainOutFile <- "Input_rain_Original.csv"
# 
# # source variables
# source_exe_path <- "C:/Program Files/eWater/Source 5.0.0.10962/"
# source_exe <- "RiverSystem.CommandLine.exe"
# project_path <- paste0(boxDir, "2021 CRAFT Barossa project/SWModel/BarossaModelFiles/Scenarios/Barossa_Sce2_stressTest1_V5.0.rsproj")
# scenario <- "Scenario 1"
# runconfiguration <-"Single analysis"
# output_file <- paste0(boxDir, "2021 CRAFT Barossa project/SWModel/BarossaModelFiles/Outputs/BaselineRun_Aug17/source_output") #changed in runSystemModelFunction
# runstarttime <- "1/1/1976"
# runendtime <- "31/12/2005"
# 
# #can create new functions for other types such as max - specify in systemArgs$outSummaryFunc below
# meanAnnual <- function(x, date) {
#   annAgg <- aggregate(x, by = list(as.POSIXlt(date)$year), FUN=sum)[["x"]]  # annual totals
#   annAggMean <- mean(annAgg)
#   return(annAggMean)
# }
# 
# systemArgs <- list()
# systemArgs$PFile <- paste0(inputDir, rainOutFile)
# systemArgs$PETFile <- paste0(inputDir, evapOutFile)
# systemArgs$startDate <- paste0(startYr, "/01/01")
# systemArgs$endDate <- paste0(endYr, "/12/31")
# systemArgs$P_colNames <- colnames(rain_data_1900to2020)
# systemArgs$PET_colNames <- c("Date", "PET")
# systemArgs$source_exe_path <- source_exe_path
# systemArgs$source_exe <- source_exe
# systemArgs$project_path <- project_path
# systemArgs$scenario <- scenario
# systemArgs$runconfiguration <- runconfiguration
# systemArgs$runstarttime <- runstarttime
# systemArgs$runendtime <- runendtime
# systemArgs$outFilename <- output_file
# # specify chosen Source results based on column named when imported to R
# systemArgs$outVarname <- c("Downstream.Flow", "Supplied.Demand.Volume", "Actual.Demand.Volume", "Actual.Evaporation.Volume", "Storage.Record", "Total.Water.Supplied", "Volume.Ordered")
# systemArgs$outSummaryFunc <- meanAnnual
# 
# metricNames <- c("Annual Yaldara Flow (ML)", "Annual Farm Dam Demand Supplied (ML)", "Annual Farm Dam Demand (ML)", "Annual Farm Dam Evaporation (ML)", "Annual Farm Dam Storage (ML)",
#                  "Annual Water Course Demand Supplied (ML)", "Annual Water Course Demand (ML)")
# 
# #################################
# 
# # To test source_wrapper if needed
# testOutput <- source_wrapper(data = sim$Rep1$Target1,
#                              systemArgs = systemArgs)
# 
# #################################
# 
# # allow output filename to be change to include rep and target info
# systemArgs$outFilename_tarRepStr='long'
# 
# systemPerf = runSystemModel(sim = sim,              # simulation; the perturbed time series
#                             systemModel = source_wrapper,      # the system model function
#                             systemArgs = systemArgs,           # argument to the system model function
#                             metrics = metricNames)             # selected performance metrics 
# 
# 
# #################################
# 
# #--------------------------------------------------------------
# # Step D: Plot the performance metrics
# #--------------------------------------------------------------
# 
# simSummary <- getSimSummary(sim)
# 
# p1 <- plotPerformanceSpace(systemPerf[1], simSummary)
# p2 <- plotPerformanceSpace(systemPerf[2], simSummary)
# p3 <- plotPerformanceSpace(systemPerf[3], simSummary)
# p4 <- plotPerformanceSpace(systemPerf[4], simSummary)
# 
# # Coverting the sum of farm dam storage to average
# systemPerfTemp=systemPerf[5]
# systemPerfTemp[[1]]=systemPerf[5][[1]]/365.25
# p5 <- plotPerformanceSpace(systemPerfTemp, simSummary)
# 
# p6 <- plotPerformanceSpace(systemPerf[6], simSummary)
# p7 <- plotPerformanceSpace(systemPerf[7], simSummary)
# 
# 
# p1_oat <- plotPerformanceOAT(systemPerf[1], simSummary)
# p2_oat <- plotPerformanceOAT(systemPerf[2], simSummary)
# p3_oat <- plotPerformanceOAT(systemPerf[3], simSummary)
# p4_oat <- plotPerformanceOAT(systemPerf[4], simSummary)
# p5_oat <- plotPerformanceOAT(systemPerfTemp, simSummary)
# p6_oat <- plotPerformanceOAT(systemPerf[6], simSummary)
# p7_oat <- plotPerformanceOAT(systemPerf[7], simSummary)
# 
# 
# 
