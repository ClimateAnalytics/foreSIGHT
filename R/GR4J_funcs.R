#####################################################################
# setup and calibrate GR4J
#' @import airGR
#' @export
setup_cal_GR4J = function(dates,P,PET,Qobs,plotResults=F){
  
  o = add_dummy_year(dates,P,PET)
  dates.new = o$dates; P.new = o$P; PET.new = o$PET
  
  ## preparation of InputsModel object
  InputsModel <- airGR::CreateInputsModel(FUN_MOD = airGR::RunModel_GR4J, DatesR = dates.new,
                                   Precip = P.new, PotEvap = PET.new)
  
  ## calibration period selection
  Ind_Run = (length(dates.new)-length(dates)+1):length(dates.new)
  IndPeriod_WarmUp = 1:(length(dates.new)-length(dates))
  
  ## preparation of RunOptions object
  RunOptions <- airGR::CreateRunOptions(FUN_MOD = airGR::RunModel_GR4J, InputsModel = InputsModel,
                                 IndPeriod_Run = Ind_Run,IndPeriod_WarmUp =IndPeriod_WarmUp )
  
  ## calibration criterion: preparation of the InputsCrit object
  InputsCrit <- airGR::CreateInputsCrit(FUN_CRIT = airGR::ErrorCrit_NSE, InputsModel = InputsModel,
                                 RunOptions = RunOptions, Obs = Qobs)
  
  ## preparation of CalibOptions object
  CalibOptions <- airGR::CreateCalibOptions(FUN_MOD = airGR::RunModel_GR4J, FUN_CALIB = airGR::Calibration_Michel)
  
  ## calibration
  OutputsCalib <- airGR::Calibration_Michel(InputsModel = InputsModel, RunOptions = RunOptions,
                                     InputsCrit = InputsCrit, CalibOptions = CalibOptions,
                                     FUN_MOD = airGR::RunModel_GR4J)
  
  ## simulation
  Param <- OutputsCalib$ParamFinalR
  OutputsModel <- airGR::RunModel_GR4J(InputsModel = InputsModel,
                                RunOptions = RunOptions, Param = Param)
  
  ## results preview
  if(plotResults){plot(OutputsModel, Qobs = Qobs[Ind_Run])}
  
  return(Param)
  
}

########################
# define year as global variable to avoid year appearing as an 
# undefined variable when performing devtools::check() 
utils::globalVariables("year")
#####################################################################

#' @import zoo airGR
#' @export
GR4J_wrapper = function(data,
                        systemArgs,
                        metrics){
  
  if (!is.null(data$P)){
    P = data$P
  } else {
    print('require P in data')
    stop()
  } 
  
  if (!is.null(data$PET)){
    PET = data$PET
  } else if (!is.null(systemArgs$PET)){
    PET = systemArgs$PET
  } else {
    print('require PET in data or systemArgs')
    stop()
  }  
  
  if (!is.null(systemArgs$dates)){
    dates = systemArgs$dates
  } else {
    print('require dates in systemArgs')
    stop()
  }
  
  o = add_dummy_year(dates,P,PET)
  dates.new = o$dates; P.new = o$P; PET.new = o$PET
  
  InputsModel <- airGR::CreateInputsModel(FUN_MOD = airGR::RunModel_GR4J, DatesR = dates.new,
                                   Precip = P.new, PotEvap = PET.new)
  
  ## calibration period selection
  Ind_Run = (length(dates.new)-length(systemArgs$dates)+1):length(dates.new)
  IndPeriod_WarmUp = 1:(length(dates.new)-length(dates))
  
  ## preparation of RunOptions object
  RunOptions <- airGR::CreateRunOptions(FUN_MOD = airGR::RunModel_GR4J, InputsModel = InputsModel,
                                 IndPeriod_Run = Ind_Run,IndPeriod_WarmUp=IndPeriod_WarmUp)
  
  ## simulation
  Param <- systemArgs$Param
  Qsim <- airGR::RunModel_GR4J(InputsModel = InputsModel,
                        RunOptions = RunOptions, Param = Param)$Qsim
  
  metricList = c()
  metricList['meanQ'] = mean(Qsim)
  metricList['P99'] = stats::quantile(Qsim,p=0.99)
  metricList['P25'] = stats::quantile(Qsim,p=0.25)
  
  # Sample daily data
  df <- data.frame(
    date = systemArgs$dates,Qsim = Qsim)  
 
  annual_data <- df %>%
    dplyr::mutate(year = lubridate::year(date)) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(
      annual_total = sum(Qsim, na.rm = TRUE))
  
  #Qsim.ann = annual_data$annual_total
  
  annual_data$rolling_3yr <- zoo::rollapply(
    annual_data$annual_total,
    width = 3,
    FUN = sum,
    align = "left",
    fill = NA
  )
  
  metricList['min3yr'] = min(annual_data$rolling_3yr,na.rm=T)

  return(metricList)
  
}

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

