rm(list=ls())

devtools::load_all()

load('data/tankDat_old.RData')
tank_obs = convert_climYMD_POSIXct(tank_obs)
save(tank_obs,file='data/tankDat.RData')


load('data/barossaDat_old.RData')
barossa_obs = convert_climYMD_POSIXct(barossa_obs)
save(barossa_obs,file='data/barossaDat.RData')



