rm(list=ls())

devtools::load_all()

clim = convert_climYMD_POSIXct(tank_obs)
  
######################################################################

# # Selected attributes
attPerturb <- c("P_day_all_tot_m", "Temp_day_all_avg_m")

# Sampling bounds and strategy
attPerturbType = "regGrid"
attPerturbSamp = c(2, 2)
attPerturbMin = c(0.8,-0.5)
attPerturbMax = c(1.2,0.5)



# Creating the exposure space
expSpace <- createExpSpace(attPerturb = attPerturb, 
                           attPerturbSamp = attPerturbSamp, 
                           attPerturbMin = attPerturbMin, 
                           attPerturbMax = attPerturbMax,
                           attPerturbType = attPerturbType)


sim <- generateScenarios(reference = clim,   # reference time series 
                         expSpace = expSpace,    # exposure space
                         controlFile = 'scaling')

pause
