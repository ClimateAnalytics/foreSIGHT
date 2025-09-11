load_camels = function(catchment,startYr,endYr,climPET=F){
  data = DroughtRisk::read_catchment_data(catchment=catchment,dataSource = 'CAMELS-AUS')
  year = as.integer(format(data$dates,'%Y'))
  keep = which((year>=startYr)&(year<=endYr))
  dates = data$dates[keep]
  times = as.POSIXct(dates,tz='UTC')
  P = data$P[keep]
  if (climPET){
    PET = data$PETclim[keep]
  } else {
    PET = data$PETobs[keep]
  }
  Qobs = data$Qobs[keep]
  return(list(times=times,P=P,PET=PET,Qobs=Qobs))
}
