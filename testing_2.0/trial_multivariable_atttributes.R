rm(list=ls())

devtools::load_all()

clim = convert_climYMD_POSIXct(tank_obs)

#calculateAttributes(clim,'P_day_all_tot')

#calculateAttributes(clim,attSel=c('P_day_all_tot','Temp_day_all_avg','mv.P.Temp_day_all_cor'))


mvFunc_cor = function(data.1,data.2){
  return(cor(data.1,data.2,use='pairwise.complete.obs'))
}

calculateAttributes(clim,attSel=c('mv.P.Temp_day_all_cor'))
