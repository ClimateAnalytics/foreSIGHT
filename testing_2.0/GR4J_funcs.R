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

setup_cal_GR4J = function(dates,P,PET,Qobs){
  
  o = add_dummy_year(dates,P,PET)
  dates.new = o$dates; P.new = o$P; PET.new = o$PET
  
  library(airGR)
  
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
  
  return(Param)
  
}


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

