rm(list=ls())

devtools::load_all()

######################
# generate monthly rainfall 

times = as.POSIXct(paste0(tank_obs$year,'/',tank_obs$month,'/',tank_obs$day),tz='UTC')

clim = list(times=times,P=tank_obs$P)

library(zoo)

P.zoo = zoo(clim$P,times)

P.zoo.mon = aggregate(P.zoo,zoo::as.yearmon,sum)

times.mon = as.POSIXct(time(P.zoo.mon),tz='UTC')

clim_mon = list(times=times.mon,
                P=coredata(P.zoo.mon))


####################

attPerturb = c("P_day_all_tot")

attPerturbType = "regGrid"
attPerturbSamp = c(1)
attPerturbMin = c(1.3)
attPerturbMax = c(1.3)

# create the exposure space
expSpace = createExpSpace(attPerturb = attPerturb,
                           attPerturbSamp = attPerturbSamp,
                           attPerturbMin = attPerturbMin,
                           attPerturbMax = attPerturbMax,
                           attPerturbType = attPerturbType)

####################

# sim = generateScenarios(reference = clim,
#                         expSpace = expSpace,
#                         controlFile = 'scaling')

sim_mon = generateScenarios(reference = clim_mon,
                        expSpace = expSpace,
                        controlFile = 'scaling')

