rm(list=ls())

devtools::load_all()

#clim = convert_climYMD_POSIXct(tank_obs)

clim = convert_climYMD_POSIXct(barossa_obs)


####################

# calculateAttributes(clim,'P_day_all_seasRatioMarMay')

####################

# #attPerturb = c("P_day_all_tot_m")
# 
# attPerturb = c("P_day_all_tot_m","P_day_all_seasRatioMarMay")
# 
# attPerturbType = "regGrid"
# attPerturbSamp = c(1,1)
# attPerturbMin = c(1.3,0.9)
# attPerturbMax = c(1.3,0.9)
# 
# # create the exposure space
# expSpace = createExpSpace(attPerturb = attPerturb,
#                           attPerturbSamp = attPerturbSamp,
#                           attPerturbMin = attPerturbMin,
#                           attPerturbMax = attPerturbMax,
#                           attPerturbType = attPerturbType)

####################

# sim_scaling = generateScenarios(reference = clim,
#                                    expSpace = expSpace,
#                                    controlFile = 'scaling')
# 
# 
# 
# 
# months = as.integer(format(clim$times,'%m'))
# keep = which(months %in% c(3,4,5))
# 
# P = sim_scaling$Rep1$Target1$P
# 
# s = 1
# 
# mean(P[,s])/mean(clim$P[,s])
# 
# SRobs = mean(clim$P[keep])/mean(clim$P)
# SRsim = mean(P[keep])/mean(P)
# SRsim/SRobs
# 
# clim_sim = clim
# clim_sim$P = P
# 
# calculateAttributes(clim_sim,'P_day_all_seasRatioMarMay')/calculateAttributes(clim,'P_day_all_seasRatioMarMay')
# 
# calculateAttributes(clim_sim,'P_day_all_tot')/calculateAttributes(clim,'P_day_all_tot')
# 
# clim_sim_ave = clim
# clim_sim_ave$P = apply(P,1,mean)
# 
# clim_ave = clim
# clim_ave$P = apply(clim$P,1,mean)
# 
# 
# calculateAttributes(clim_sim_ave,'P_day_all_seasRatioMarMay')/calculateAttributes(clim_ave,'P_day_all_seasRatioMarMay')
# calculateAttributes(clim_sim_ave,'P_day_all_tot')/calculateAttributes(clim_ave,'P_day_all_tot')

####################

# generate monthly rainfall 

library(zoo)

P.zoo = zoo(clim$P,clim$times)

P.zoo.mon = aggregate(P.zoo,zoo::as.yearmon,sum)

times.mon = time(P.zoo.mon)

clim_mon = list(times=as.POSIXct(times.mon,tz = 'UTC'),
                P=coredata(P.zoo.mon))

####################

attPerturb = c("P_month_all_tot_m","P_month_all_seasRatioMarMay")

attPerturbType = "regGrid"
attPerturbSamp = c(1,1)
attPerturbMin = c(1.3,0.9)
attPerturbMax = c(1.3,0.9)

# create the exposure space
expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType)



sim_mon_scaling = generateScenarios(reference = clim_mon,
                                expSpace = expSpace,
                                controlFile = 'scaling')

months = as.integer(format(clim_mon$times,'%m'))
keep = which(months %in% c(3,4,5))

P = sim_mon_scaling$Rep1$Target1$P

s = 1

mean(P[,s])/mean(clim_mon$P[,s])

SRobs = mean(clim_mon$P[keep])/mean(clim_mon$P)
SRsim = mean(P[keep])/mean(P)
SRsim/SRobs
