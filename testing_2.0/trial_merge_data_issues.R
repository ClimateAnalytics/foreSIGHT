rm(list=ls())

devtools::load_all()

#clim = convert_climYMD_POSIXct(tank_obs)

clim = convert_climYMD_POSIXct(barossa_obs)


####################

calculateAttributes(clim,'P_day_Jan_tot')

