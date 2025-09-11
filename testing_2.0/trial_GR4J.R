rm(list=ls())

#####################################################################
# setup data

devtools::load_all('/Users/a1065639/Work/DroughtRisk/')

#catchment = '410057'
catchment = 'A5030502'

#startYr = 1976
#endYr = 2005

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

add_dummy_year = function(dates,P,PET,warmupYrs=1){
  year1 = as.integer(format(dates,'%Y'))[1]
  year0 = year1-warmupYrs
  dates.tmp = seq(as.POSIXct(paste0(year0,'/01/01'),tz='UTC'),
                  as.POSIXct(paste0(year0,'/12/31'),tz='UTC'),
                  by='days')
  P.tmp = P[1:length(dates.tmp)]
  PET.tmp = PET[1:length(dates.tmp)]
  
  dates.new = c(dates.tmp,dates)
  P.new = c(P.tmp,P)
  PET.new = c(PET.tmp,PET)
  
  return(list(dates=dates.new,
              P = P.new,
              PET = PET.new))
  
}

#####################################################################
# setup and calibrate GR4J

o = add_dummy_year(dates,P,PET)
dates.new = o$dates; P.new = o$P; PET.new = o$PET

library(airGR)

## loading catchment data
#data(L0123001)

## preparation of InputsModel object
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = dates.new,
                                 Precip = P.new, PotEvap = PET.new)

## calibration period selection
# Ind_Run <- seq(which(format(dates.new, format = "%Y-%m-%d")==paste0(startYr,'-01-01')),
#                which(format(dates.new, format = "%Y-%m-%d")==paste0(endYr,'-12-31')))
Ind_Run = (length(dates.new)-length(dates)+1):length(dates.new)
IndPeriod_WarmUp = 1:(length(dates.new)-length(dates))

## preparation of RunOptions object
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J, InputsModel = InputsModel,
                               IndPeriod_Run = Ind_Run,IndPeriod_WarmUp =IndPeriod_WarmUp )

## calibration criterion: preparation of the InputsCrit object
InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModel,
                               RunOptions = RunOptions, Obs = Qobs)

## preparation of CalibOptions object
CalibOptions <- CreateCalibOptions(FUN_MOD = RunModel_GR4J, FUN_CALIB = Calibration_Michel)

## calibration
OutputsCalib <- Calibration_Michel(InputsModel = InputsModel, RunOptions = RunOptions,
                                   InputsCrit = InputsCrit, CalibOptions = CalibOptions,
                                   FUN_MOD = RunModel_GR4J)

## simulation
Param <- OutputsCalib$ParamFinalR
OutputsModel <- RunModel_GR4J(InputsModel = InputsModel,
                              RunOptions = RunOptions, Param = Param)

## results preview
plot(OutputsModel, Qobs = data$Qobs[Ind_Run])


#####################################################################

GR4J_wrapper = function(data,
                        systemArgs,
                        metrics){
  
  o = add_dummy_year(systemArgs$dates,data$P,data$PET)
  dates.new = o$dates; P.new = o$P; PET.new = o$PET
  
  InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = dates.new,
                                   Precip = P.new, PotEvap = PET.new)
  
  ## calibration period selection
  Ind_Run = (length(dates.new)-length(systemArgs$dates)+1):length(dates.new)
  IndPeriod_WarmUp = 1:(length(dates.new)-length(dates))
  
  ## preparation of RunOptions object
  RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J, InputsModel = InputsModel,
                                 IndPeriod_Run = Ind_Run,IndPeriod_WarmUp=IndPeriod_WarmUp)
  
  ## simulation
  Param <- systemArgs$Param
  Qsim <- RunModel_GR4J(InputsModel = InputsModel,
                        RunOptions = RunOptions, Param = Param)$Qsim
  
  metricList = c()
  metricList['meanQ'] = mean(Qsim)
  metricList['P99'] = quantile(Qsim,p=0.99)
  metricList['P25'] = quantile(Qsim,p=0.25)
  
  return(metricList)
  
}

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

attPerturb = c('P_day_all_tot_m',"PET_day_all_avg")

attPerturbType = "regGrid"
attPerturbSamp = c(5,5)
attPerturbMin = c(0.8,1)
attPerturbMax = c(1.2,1.4)

expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType)

sim_scaling = generateScenarios(reference = clim_ref,
                                expSpace = expSpace,
                                controlFile = 'scaling')

systemArgs = list(dates=dates,Param=Param)

sysOut = runSystemModel(sim=sim_scaling,systemModel = GR4J_wrapper,systemArgs = systemArgs,metrics = c('meanQ','P99','P25'))

plotPerformanceSpace(performance = sysOut,sim=sim_scaling)

plotPerformanceOAT(performance = sysOut,sim=sim_scaling,metric = 'P25')
plotPerformanceOAT(performance = sysOut,sim=sim_scaling,metric = 'P99')

#####################################################################

# attPerturb = c('P_day_all_tot_m',"P_day_all_seasRatioJunAug","PET_day_all_avg")
# #attPerturb = c("P_day_all_seasRatioJunAug","PET_day_all_avg")
# #attPerturb = c('P_day_all_tot_m',"PET_day_all_avg")
# 
# attPerturbType = "regGrid"
# attPerturbSamp = c(5,5,1)
# attPerturbMin = c(0.7,1.1,1)
# attPerturbMax = c(0.8,1.2,1)

attPerturb <- c("PET_day_all_avg_m", "P_day_all_tot_m","P_day_all_seasRatioMarMay") # note that updates to attribute manager allow customisation of months used for seasRatio 

# specify perturbation type and minimum-maximum ranges of the perturbed attributes
attPerturbType <- "regGrid"
#attPerturbSamp <- c(3,5,5)     
#attPerturbMin <- c(1, 0.7,0.7)
#attPerturbMax <- c(1.2, 1.3,1.3)
attPerturbSamp <- c(1,7,7)     
attPerturbMin <- c(1, 0.8,0.5)
attPerturbMax <- c(1, 1.2,1.8)

expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType)

sim_seasScaling = generateScenarios(reference = clim_ref,
                                    expSpace = expSpace,
                                    controlFile = 'scaling')

systemArgs = list(dates=dates,Param=Param)

sysOut = runSystemModel(sim=sim_seasScaling,systemModel = GR4J_wrapper,systemArgs = systemArgs,metrics = 'meanQ')

plotPerformanceSpace(performance = sysOut,sim=sim_seasScaling,attX = "P_day_all_tot_m",attY="P_day_all_seasRatioMarMay")

######################################################################

func_cv = function(data){
  m = mean(data,na.rm=T)
  if (m==0){
    cv = 9999.
  } else {
    cv = sd(data,na.rm=T)/m
  }
  return(cv)
}

mvFunc_cvWetDay = function(data.1,data.2){
  return(func_cv(data.1[data.2>0]))
}

mvFunc_cvDryDay = function(data.1,data.2){
  return(func_cv(data.1[data.2==0]))
}

mvFunc_avgWetDay = function(data.1,data.2){
  return(mean(data.1[data.2>0],na.rm=T))
}

mvFunc_avgDryDay = function(data.1,data.2){
  return(mean(data.1[data.2==0],na.rm=T))
}

######################################################################

modelSelection = list()

modelSelection$modelType = list()
modelSelection$modelType$P = "latent"
modelSelection$modelType$PET = "wgenLM"

modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "seas"
modelSelection$modelParameterVariation$PET = "harWD"

modelSelection[["optimisationArguments"]] <- list()
modelSelection[["optimisationArguments"]][["nMultiStart"]] <- 5
modelSelection[["optimisationArguments"]][['RGN.control']]=list(iterMax=50)

modelSelection[["penaltyAttributes"]] <- c("P_day_all_tot", "P_day_all_xP99overPave",
                                           "P_day_all_nWet", "P_day_all_avgDSD",
                                           "PET_day_all_avg","PET_day_all_cor","PET_day_all_cv")
modelSelection[["penaltyWeights"]] = rep(2,7)


modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")
write(modelSelectionJSON, file = controlFile)

####################

attPerturb = c('P_day_all_tot')
attHold = c('P_day_all_avgDSD','P_day_all_xP99overPave','P_day_all_nWet',
            "PET_day_all_avg","PET_day_all_cor","PET_day_all_cv")

attPerturbType = "regGrid"
attPerturbSamp = c(5)
attPerturbMin = c(0.8)
attPerturbMax = c(1.2)

#########################

targetTypes = list()

# create the exposure space
expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType,
                          attHold = attHold,
                          targetTypes = targetTypes)

attTied = list(wDdD=c("PET_day_all_avg","PET_day_all_cv"),
               seas='allTargets')

expSpace = tieAttributes(expSpace,attTied)

# i=which(grepl(pattern = 'xP99overPave',x = colnames(expSpace$targetMat)))
# expSpace$targetMat[,i] = 1/expSpace$targetMat[,'P_day_all_tot']


time.1 = Sys.time()
sim_stoch = generateScenarios(reference = clim_ref,
                              expSpace = expSpace,
                              controlFile = controlFile,
                              seedID = 1,
                              numReplicates = 5)
time.2 = Sys.time()
print(time.2-time.1)

plotScenarios(sim_stoch)

sysOutSim = runSystemModel(sim=sim_stoch,systemModel = GR4J_wrapper,systemArgs = systemArgs,metrics = c('meanQ','P99','P25'))

sysOutClim = GR4J_wrapper(data = clim_ref, systemArgs = systemArgs,metrics = c('meanQ','P99','P25'))

plotPerformanceOAT(performance = sysOutSim, sim=sim_stoch, metric = 'meanQ')
plotPerformanceOAT(performance = sysOutSim, sim=sim_stoch, metric = 'P99')
plotPerformanceOAT(performance = sysOutSim, sim=sim_stoch, metric = 'P25')

pause

#########################





modelSelectionA = list()

modelSelectionA$modelType = list()
modelSelectionA[["optimisationArguments"]] <- list()
modelSelectionA[["optimisationArguments"]][["nMultiStart"]] <- 5
modelSelectionA[["optimisationArguments"]][['RGN.control']]=list(iterMax=50)

modelSelectionJSONA = jsonlite::toJSON(modelSelectionA, pretty = TRUE, auto_unbox = TRUE)
controlFileA = paste0(tempdir(), "\\eg_controlFileA.json")
write(modelSelectionJSONA, file = controlFileA)

attPerturb <- c("P_day_all_tot_m")
attHold <- c("P_day_all_seasRatio", "P_day_all_nWet_m", "P_day_all_maxDSD_m", "P_day_all_maxWSD_m")
attPerturbType = "regGrid"
attPerturbSamp = c(5)
attPerturbMin = c(0.8)
attPerturbMax = c(1.2)
expSpace <- createExpSpace(attPerturb = attPerturb, 
                           attPerturbSamp = attPerturbSamp, 
                           attPerturbMin = attPerturbMin, 
                           attPerturbMax = attPerturbMax,
                           attPerturbType = attPerturbType, 
                           attHold = attHold)
data(tankDat)        # reference data
simHold <- generateScenarios(reference = clim_ref, 
                             expSpace = expSpace,
                             controlFile = controlFileA,
                             seedID = 1,
                             numReplicates = 5) # simulation

plotScenarios(simHold)

sysOutSimHold = runSystemModel(sim=simHold,systemModel = GR4J_wrapper,systemArgs = systemArgs,metrics = c('meanQ','P99','P25'))

sysOutClim = GR4J_wrapper(data = clim_ref, systemArgs = systemArgs,metrics = c('meanQ','P99','P25'))

plotPerformanceOAT(performance = sysOutSimHold, sim=sysOutSimHold, metric = 'meanQ')
plotPerformanceOAT(performance = sysOutSimHold, sim=sysOutSimHold, metric = 'P99')
plotPerformanceOAT(performance = sysOutSimHold, sim=sysOutSimHold, metric = 'P25')

attSel = colnames(simHold$expSpace$targetMat)
par(mfrow=c(4,4))
P = calcPerformanceAttributes(clim=clim_ref,sim=simHold,attSel=attSel)
for (att in names(P)){
  plotPerformanceOAT(P, simHold, metric=att,plotType='base',attSel=attPerturb)
}

attSel = colnames(sim_stoch$expSpace$targetMat)
attSel = attSel[startsWith(x = attSel,prefix='P_')]
par(mfrow=c(4,4))
P = calcPerformanceAttributes(clim=clim_ref,sim=simHold,attSel=attSel)
for (att in names(P)){
  plotPerformanceOAT(P, simHold, metric=att,plotType='base',attSel=attPerturb)
}

#########################


attPerturb = c("PET_day_all_avg")
attHold = c('P_day_all_tot_m','P_day_all_avgDSD','P_day_all_xP99overPave','P_day_all_nWet',
            "PET_day_all_cor","PET_day_all_cv")

attPerturbType = "regGrid"
#attPerturbSamp = c(5)
#attPerturbMin = c(0.7)
#attPerturbMax = c(1.1)

attPerturbSamp = c(1)
attPerturbMin = c(1)
attPerturbMax = c(1)

#########################

# attPerturbSamp = c(1)
# attPerturbMin = c(1)
# attPerturbMax = c(1)

# create the exposure space
expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType,
                          attHold = attHold)#,
#                          attTied = c(attPerturb,attHold),
#                          tieType = 'seas')

attTied = list(wDdD=c("PET_day_all_avg","PET_day_all_cv"),
               seas='allTargets')

expSpace = tieAttributes(expSpace,attTied)


#attTied = c("PET_day_all_avg","PET_day_all_xP90")
#tieType = 'wDdD'

# attPerturb = expSpace$attPerturb
# attHold = expSpace$attHold
# 
# attTied = c(attPerturb,attHold)

# attTied = c("PET_day_all_avg","PET_day_all_xP90")
# 
# if (!all(attTied%in%c(attPerturb,attHold))){stop("must have tied attributes in perturb/hold atts")}
# 
# if (!is.null(attTied)){
#   if(tieType=='wDdD'){
#     for (att in attTied){
#       i=which(colnames(expSpace$targetMat)==att)
#       var = get.attribute.varType(att)
#       if(grepl('/',var)){stop("can't have multivariable tied attributes in perturb/hold atts")}
#       var.new = paste0(var,'.P')
#       for (cond in c('WetDay','DryDay')){
#         att.cond = gsub(var,var.new,att)
#         att.cond = paste0('mv.',att.cond,cond)
#         if (att.cond%in%expSpace$attTied){
#           expSpace$targetMat[att.cond] = expSpace$targetMat[att.cond]*expSpace$targetMat[att]
#         } else {
#           expSpace$targetMat[att.cond] = expSpace$targetMat[att]
#           expSpace$attTied = c(expSpace$attTied,att.cond)
#           expSpace$targetType = c(expSpace$targetType,expSpace$targetType[i])
#         }
#       }
#     }
#   }
# }
# 
# #################
# 
# #attTied = c(attPerturb,attHold)
# tieType = 'seas'
# 
# #attPerturb = expSpace$attPerturb
# #attHold = expSpace$attHold
# attTied = colnames(expSpace$targetMat)
# #if (!all(attTied%in%c(attPerturb,attHold))){stop("must have tied attributes in perturb/hold atts")}
# 
# if (!is.null(attTied)){
#   if(tieType=='seas'){
#     for (att in attTied){
#       i=which(colnames(expSpace$targetMat)==att)
#       for (seas in c('DJF','MAM','JJA','SON')){
#         att.seas = gsub('all',seas,att)
#         if (att.seas%in%expSpace$attTied){
#           expSpace$targetMat[att.seas] = expSpace$targetMat[att.seas]*expSpace$targetMat[att]
#         } else {
#           expSpace$targetMat[att.seas] = expSpace$targetMat[att]
#           expSpace$attTied = c(expSpace$attTied,att.seas)
#           expSpace$targetType = c(expSpace$targetType,expSpace$targetType[i])
#         }
#       }
#     }
#   }
# }

####################

# calculateAttributes(clim_ref,attSel = colnames(expSpace$targetMat))
# 
# P = clim_ref$P
# PET = clim_ref$PET
# month = as.integer(format(clim_ref$times,'%m'))
# keep = month%in%c(9,10,11)
# P[!keep] = NA
# PET[!keep] = NA
# 
# func_cor(PET)
# 
# pause

time.1 = Sys.time()
sim_stoch = generateScenarios(reference = clim_ref,
                              expSpace = expSpace,
                              controlFile = controlFile,
                              seedID = 1,
                              numReplicates = 1)
time.2 = Sys.time()
print(time.2-time.1)

plotScenarios(sim_stoch)

attSel = colnames(expSpace$targetMat)

P = calcPerformanceAttributes(clim=clim_ref,sim=sim_stoch,attSel=attSel)

par(mfrow=c(5,5),mar=c(4,5,2,1))
for (att in names(P)){
  # plotPerformanceSpace(P, sim, metric=metric)
  plotPerformanceOAT(P, sim_stoch, metric=att,plotType='base',attSel=attPerturb)
}

####################

