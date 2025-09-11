rm(list=ls())

########################

foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT/'
droughRiskDir = 'C:/Users/a1065639/Work/DroughtRisk/'
dataDir.Loxton = 'C:/Users/a1065639/Box/2025_DEW_foreSIGHT/Data/'
dataDir.Barossa = 'C:/Users/a1065639/Box/2021 CRAFT Barossa project/External datasets/Hydroclimate (instrumental)/'

# foreSIGHTDir = '../'
# droughRiskDir = '/hpcfs/users/a1065639/git/DroughtRisk_Feb_2024_paper_shared/DroughtRisk/'
# dataDir.Loxton = '../../Data/'
# dataDir.Barossa = '../../Data/'

devtools::load_all(foreSIGHTDir)
devtools::load_all(droughRiskDir)

runDirname = paste0(foreSIGHTDir,'testing_2.0/')

########################

# location
#catchment = 'Loxton'
#catchment = '410057' # Lacmalac
catchment = 'A5030502' # Scott Creek

#catchment = '922101B' # Coen River at Racecourse 
# catchment = '116014A' # Wild River at Silver Valley 
# catchment = '419005' # Namoi River at North Cuerindi 
#catchment = '412066' # Abercrombie River at Hadley No. 2 
# catchment = '410734' # Queanbeyan River at Tinderry 
# catchment = '410730' # Cotter River at Gingera  
# catchment = '401012' # Murray River at Biggara 
# catchment = '61' # Hellyer River at Guilford Junction 
# catchment = '231213' # Lerderderg River at O’Brien Crossing 
# catchment = 'A5130501' # Rocky River upstream Gorge Falls 
# catchment = '606001' # Deep River at Teds Pool 

# catchment = 'Barossa.1'

# period
startYr = 1976
endYr = 2005
# startYr = 1980
# endYr = 1989

# perturberd attribute
#attPerturb = "P_day_all_tot"
#attPerturb = "P_day_all_xP99overPave"
#attPerturb = "P_day_all_avgDSD"
#attPerturb = "PET_day_all_avg"
attPerturb = "P_day_all_seasRatioMarMay"

numReplicates = 3
cores = 1

# args = commandArgs(trailingOnly=TRUE)
# catchment = args[1]
# attPerturb = args[2]
# startYr = args[3]
# endYr = args[4]
# numReplicates = args[5]
# cores = args[6]

########################

if (catchment=='Loxton'){
  fname = paste0(dataDir.Loxton,'24007_LOXTON_Dec2023.csv')
  d = read.csv(fname,header=T,skip=4017) # note there is different format for dates pre and post 1900. have chosen to not read in pre-1900 data.
  date = as.Date(d[,1],format = '%d/%m/%Y')
  P = d[,2]
  PET = d[,3]
  # select time periods of data to use for simulations
  year = as.integer(format(date,'%Y'))
  keep = which((year>=startYr)&(year<=endYr))
  times = as.POSIXct(date[keep],tz='UTC')
  P = P[keep]
  PET = PET[keep]
  Qobs = NULL
} else if (catchment=='Barossa.1'){
site_keep = 1:1
	#site_keep = 1:14
  evapFile <- "evap_1900to2020.Rdata"
  rainFile <- "rain_data_1900to2020_heggiesInfilled.Rdata"
  load(paste0(dataDir.Barossa, evapFile))
  evap.date = as.Date(evap_1900to2020[,1],format = '%d/%m/%Y')
  evap.data = evap_1900to2020[,2]
  load(paste0(dataDir.Barossa, rainFile))
  nRain = dim(rain_data_1900to2020)[2]-1
  rain.date = as.Date(rain_data_1900to2020[,1],format = '%d/%m/%Y')
  rain.data = rain_data_1900to2020[,2:(nRain+1)]
  if (!any(rain.date!=evap.date)){
    date=rain.date
  } else {
    stop('dates for P and PET file differ')
  }
  year = as.integer(format(date,'%Y'))
  keep = which((year>=startYr)&(year<=endYr))
  times = as.POSIXct(date[keep],tz='UTC')
  P = rain.data[keep,site_keep]
  PET = evap.data[keep]
  Qobs = NULL
} else {
  data = DroughtRisk::read_catchment_data(catchment=catchment,dataSource = 'CAMELS-AUS')
  year = as.integer(format(data$dates,'%Y'))
  keep = which((year>=startYr)&(year<=endYr))
  dates = data$dates[keep]
  times = as.POSIXct(dates,tz='UTC')
  P = data$P[keep]
  PET = data$PETobs[keep]
  Qobs = data$Qobs[keep]
}

# create reference data 
clim_ref <- list(times = times,
                 P = P,
                 PET = PET)

########################

mvFunc_cor = function(data.1,data.2){
  return(cor(data.1,data.2,use='pairwise.complete.obs'))
}

mvFunc_avgWetDay = function(data.1,data.2){
  return(mean(data.1[data.2>0],na.rm=T))
}

mvFunc_avgDryDay = function(data.1,data.2){
  return(mean(data.1[data.2==0],na.rm=T))
}

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

######################################################################

modelSelection = list()

modelSelection$modelType = list()
modelSelection$modelType$P = "latent"
modelSelection$modelType$PET = "wgenLM"

modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "seas"
modelSelection$modelParameterVariation$PET = "harWD"

modelSelection$postProcessing = list(P=list())
modelSelection$postProcessing$P$types = c('annVar','scaleExtremesSeas')

modelSelection[["optimisationArguments"]] <- list()
modelSelection[["optimisationArguments"]][["nMultiStart"]] <- 5
modelSelection[["optimisationArguments"]][["OFtol"]] <- 0.02
modelSelection[["optimisationArguments"]][['RGN.control']]=list(iterMax=50)

#modelSelection[["penaltyAttributes"]] <- c("P_day_all_tot", "P_day_all_xP99overPave",
#                                           "P_day_all_nWet", "P_day_all_avgDSD",
#                                           "P_day_all_seasRatioMarMay","P_day_all_tot_cv",
#                                           "PET_day_all_avg","PET_day_all_cor","PET_day_all_cv")
modelSelection[["penaltyAttributes"]] <- c("P_day_all_tot", "P_day_all_P99",
                                           "P_day_all_nWet", "P_day_all_avgDSD",
                                           "P_day_all_seasRatioMarMay","P_day_all_tot_cv",
                                           "PET_day_all_avg","PET_day_all_cor","PET_day_all_cv")
modelSelection[["penaltyWeights"]] = rep(3,length(modelSelection[["penaltyAttributes"]]))

modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")
write(modelSelectionJSON, file = controlFile)

####################

attPerturbType = "OAT"
#attPerturb = c('P_day_all_tot','P_day_all_xP99overPave','P_day_all_avgDSD')
# attPerturb = c('P_day_all_tot','P_day_all_P99','P_day_all_avgDSD')
# attPerturbSamp = c(5,5,5)
# attPerturbMin = c(0.8,0.9,1)
# attPerturbMax = c(1.2,1.3,1.4)
# attHold = c('P_day_all_nWet',"P_day_all_seasRatioMarMay",'P_day_all_tot_cv',
#             "PET_day_all_avg","PET_day_all_cor","PET_day_all_cv")

attPerturb = c('P_day_all_avgDSD')
attPerturbSamp = c(5)
attPerturbMin = c(1)
attPerturbMax = c(1.4)
attHold = c('P_day_all_tot','P_day_all_P99','P_day_all_nWet',"P_day_all_seasRatioMarMay",'P_day_all_tot_cv',
            "PET_day_all_avg","PET_day_all_cor","PET_day_all_cv")


attsAll = c(attPerturb,attHold)

expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType,
                          attHold = attHold)


attTied = list(wDdD=c("PET_day_all_avg","PET_day_all_cv"),
               seas=attsAll[!attsAll%in%c('P_day_all_seasRatioMarMay','P_day_all_tot_cv')])

# if (attPerturb=='P_day_all_seasRatioMarMay'){
#   exclude.seas=c('P_day_MAM_tot','P_day_SON_tot')
# } else {
  exclude.seas=NULL
# }
expSpace = tieAttributes(expSpace=expSpace,attTied=attTied,exclude.seas=exclude.seas)

# #if (attPerturb=='P_day_all_tot'){
#   i=which(grepl(pattern = 'xP99overPave',x = colnames(expSpace$targetMat)))
#   expSpace$targetMat[,i] = 1/expSpace$targetMat[,'P_day_all_tot']
# #}

time.1 = Sys.time()
sim_stoch = generateScenarios(reference = clim_ref,
                              expSpace = expSpace,
                              controlFile = controlFile,
                              seedID = 1,
                              numReplicates = numReplicates,
                              cores = cores)
time.2 = Sys.time()
print(time.2-time.1)

pause
####################

attsAll = c('P_day_all_avgDSD','P_day_all_tot','P_day_all_xP99overPave',
            'P_day_all_nWet',"P_day_all_seasRatioMarMay",'P_day_all_tot_cv',
            "PET_day_all_avg","PET_day_all_cor","PET_day_all_cv")

####################

attPerturbType = "regGrid"
attPerturbSamp = c(5)
attHold = attsAll[attsAll!=attPerturb]
if (attPerturb=='P_day_all_avgDSD'){
  attPerturbMin = c(1.)
  attPerturbMax = c(1.4)
} else if (attPerturb=="PET_day_all_avg"){
  attPerturbMin = c(1)
  attPerturbMax = c(1.2)
} else if (attPerturb=='P_day_all_tot'){
  attPerturbMin = c(0.8)
  attPerturbMax = c(1.2)
} else if (attPerturb=='P_day_all_xP99overPave'){
  attPerturbMin = c(0.9)
  attPerturbMax = c(1.3)
} else if (attPerturb=='P_day_all_seasRatioMarMay'){
  attPerturbMin = c(0.5)
  attPerturbMax = c(1.5)
}

print('not perturbing climate')
attPerturbSamp = c(1)
attPerturbMin = c(1.)
attPerturbMax = c(1)
  
##########################

# create the exposure space
expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType,
                          attHold = attHold)

attTied = list(wDdD=c("PET_day_all_avg","PET_day_all_cv"),
               seas=attsAll[!attsAll%in%c('P_day_all_seasRatioMarMay','P_day_all_tot_cv')])

if (attPerturb=='P_day_all_seasRatioMarMay'){
  exclude.seas=c('P_day_MAM_tot','P_day_SON_tot')
} else {
  exclude.seas=NULL
}
expSpace = tieAttributes(expSpace=expSpace,attTied=attTied,exclude.seas=exclude.seas)

if (attPerturb=='P_day_all_tot'){
  i=which(grepl(pattern = 'xP99overPave',x = colnames(expSpace$targetMat)))
  expSpace$targetMat[,i] = 1/expSpace$targetMat[,'P_day_all_tot']
}



# f = 1.3
# 
# o = calculateAttributes(clim_ref,c('P_day_all_tot','P_day_JJA_tot'))
# 
# g = 1+(1-f)*o['P_day_JJA_tot']/(o['P_day_all_tot']-o['P_day_JJA_tot'])
# 
# i = which(colnames(expSpace$targetMat)=="P_day_JJA_tot" )
# 
# expSpace$targetMat[,i] = expSpace$targetMat[,i]*f
# 
# i = which(colnames(expSpace$targetMat)%in%c("P_day_DJF_tot","P_day_MAM_tot","P_day_SON_tot"))
# expSpace$targetMat[,i] = expSpace$targetMat[,i]*g



time.1 = Sys.time()
sim_stoch = generateScenarios(reference = clim_ref,
                              expSpace = expSpace,
                              controlFile = controlFile,
                              seedID = 1,
                              numReplicates = numReplicates,
                              cores = cores)
time.2 = Sys.time()
print(time.2-time.1)

pause

sim_stoch_shuffle = shuffle_sim(sim=sim_stoch,clim=clim_ref,attPerturb='P_day_all_tot_dwellTime',
             targetVals=c(1,1.25,1.5))

# ####################################
# 
# target_dwellTime_shuffle = function(times,
#                                     P.sim,
#                                     target.dwellTime,
#                                     seed=1,
#                                     annAR1coeffList = seq(-0.2,0.9,0.01)){
# 
#   clim_sim = list(times=times,P=P.sim)
#   dwellTime.sim = calculateAttributes(clim_sim,'P_day_all_tot_dwellTime')
#   
#   dwellTime.sim.multi = c()
#   for (i in 1:length(annAR1coeffList)){
#     annAR1coeff = annAR1coeffList[i]
#     P.new = pp.annShuffle(P.sim,times,annAR1coeff,seed)
#     clim_sim = list(times=times,P=P.new)
#     dwellTime.sim.multi[i] = calculateAttributes(clim_sim,'P_day_all_tot_dwellTime')
#   }
#  
#   abs_diff = abs(dwellTime.sim.multi - dwellTime.target)
#   i = floor(median(which(abs_diff == min(abs_diff))))
#   
#   annAR1coeff = annAR1coeffList[i]
#   P.targ = pp.annShuffle(P,times,annAR1coeff,seed)
#   clim_sim_targ = list(times=times,P=P.targ)
#   
#   browser()
#   
#   return(P.targ)
#   
# }
# 
# dwellTime.clim = calculateAttributes(clim_ref,'P_day_all_tot_dwellTime')
# dwellTime.target = dwellTime.clim*1.5
# 
# P.targ = target_dwellTime_shuffle(times=sim_stoch$simDates,
#                                     P.sim=sim_stoch$Rep1$Target1$P$sim,
#                                     target.dwellTime=dwellTime.target)
# 
# # pause
# # 
# # min5yr.clim = calculateAttributes(clim_ref,'P_day_all_tot_min5yr')
# # max5yr.clim = calculateAttributes(clim_ref,'P_day_all_tot_max5yr')
# # dwellTime.clim = calculateAttributes(clim_ref,'P_day_all_tot_dwellTime')
# # 
# # P = sim_stoch$Rep1$Target1$P$sim
# # times = sim_stoch$simDates
# # clim_sim = list(times=times,P=P)
# # min5yr.sim = calculateAttributes(clim_sim,'P_day_all_tot_min5yr')
# # max5yr.sim = calculateAttributes(clim_sim,'P_day_all_tot_max5yr')
# # dwellTime.sim = calculateAttributes(clim_sim,'P_day_all_tot_dwellTime')
# # 
# # seed = 1
# # iyy = NULL
# # 
# # # annAR1coeff = 0.
# # # P.0 = pp.annShuffle(P,times,annAR1coeff,seed,iyy)
# # # clim_sim = list(times=times,P=P.0)
# # # calculateAttributes(clim_sim,'P_day_all_tot_min5yr')
# # # calculateAttributes(clim_sim,'P_day_all_tot_max5yr')
# # # 
# # # 
# # # annAR1coeff = 0.8
# # # P.0.8 = pp.annShuffle(P,times,annAR1coeff,seed,iyy)
# # # clim_sim = list(times=times,P=P.0.8)
# # # calculateAttributes(clim_sim,'P_day_all_tot_min5yr')
# # # calculateAttributes(clim_sim,'P_day_all_tot_max5yr')
# # 
# # 
# # annAR1coeffList = seq(-0.2,0.9,0.01)
# # min5yr.sim.multi = max5yr.sim.multi = dwellTime.sim.multi = c()
# # for (i in 1:length(annAR1coeffList)){
# #   annAR1coeff = annAR1coeffList[i]
# #   P.new = pp.annShuffle(P,times,annAR1coeff,seed,iyy)
# #   clim_sim = list(times=times,P=P.new)
# #   min5yr.sim.multi[i] = calculateAttributes(clim_sim,'P_day_all_tot_min5yr')
# #   max5yr.sim.multi[i] = calculateAttributes(clim_sim,'P_day_all_tot_max5yr')
# #   dwellTime.sim.multi[i] = calculateAttributes(clim_sim,'P_day_all_tot_dwellTime')
# # }
# # 
# # par(mfrow=c(1,1))
# # ylim = c(min(min5yr.sim.multi,max5yr.sim.multi),max(min5yr.sim.multi,max5yr.sim.multi))
# # plot(annAR1coeffList,min5yr.sim.multi)
# # #lines(annAR1coeffList,max5yr.sim.multi)
# # abline(h=min5yr.clim,col='blue')
# # abline(h=min5yr.sim,col='red')
# # abline(h=0.8*min5yr.clim,col='green')
# # 
# # dwellTime.target = dwellTime.clim*1.5
# # abs_diff = abs(dwellTime.sim.multi - dwellTime.target)
# # i = floor(median(which(abs_diff == min(abs_diff))))
# # 
# # # min5yr.target = min5yr.clim*0.8
# # # abs_diff = abs(min5yr.sim.multi - min5yr.target)
# # # i = which(abs_diff == min(abs_diff))
# # 
# # annAR1coeff = annAR1coeffList[i]
# # P.targ = pp.annShuffle(P,times,annAR1coeff,seed,iyy)
# # clim_sim_targ = list(times=times,P=P.targ)
# # 
# # 
# # par(mfrow=c(1,1))
# # ylim = c(min(dwellTime.sim.multi,dwellTime.clim,dwellTime.sim),
# #          max(dwellTime.sim.multi,dwellTime.clim,dwellTime.sim))
# # plot(annAR1coeffList,dwellTime.sim.multi,ylim=ylim)
# # #lines(annAR1coeffList,max5yr.sim.multi)
# # abline(h=dwellTime.clim,col='blue')
# # abline(h=dwellTime.sim,col='red')
# # abline(h=1.5*dwellTime.clim,col='green')
# # 
# # 
# # 
# # # set.seed(1)
# # # dates <- seq(as.Date("2010-01-01"), as.Date("2020-12-31"), by = "day")
# # # values <- rnorm(length(dates))
# 
# df.clim = data.frame(date = clim_ref$times, value = clim_ref$P)
# df.clim$year <- format(df.clim$date, "%Y")
# annual.clim <- aggregate(value ~ year, data = df.clim, FUN = mean)
# 
# df.sim = data.frame(date = sim_stoch$simDates, value = sim_stoch$Rep1$Target1$P$sim)
# df.sim$year <- format(df.sim$date, "%Y")
# annual.sim <- aggregate(value ~ year, data = df.sim, FUN = mean)
# 
# df.sim.targ = data.frame(date = sim_stoch$simDates, value = P.targ)
# df.sim.targ$year <- format(df.sim$date, "%Y")
# annual.sim.targ <- aggregate(value ~ year, data = df.sim.targ, FUN = mean)
# 
# 
# plot(annual.clim$year,annual.clim$value,type='o')
# abline(h=median(annual.clim$value),lty=2)
# lines(annual.sim$year,annual.sim$value,col='red',type='o')
# abline(h=median(annual.sim$value),col='red',lty=2)
# lines(annual.sim.targ$year,annual.sim.targ$value,col='blue',type='o')
# abline(h=median(annual.sim.targ$value),col='blue',lty=2)

####################################



fname = paste0(runDirname,catchment,'_',attPerturb,'.RData') 
save.image(file=fname)

fig.fname = paste0(runDirname,catchment,'_',attPerturb,'.pdf') 
pdf(fig.fname)

plotScenarios(sim_stoch)

attSel = colnames(expSpace$targetMat)
attSel = c(attSel,'P_day_all_tot_dwellTime')

P = calcPerformanceAttributes(clim=clim_ref,sim=sim_stoch,attSel=attSel)


par(mfrow=c(5,5),mar=c(4,5,2,1))
for (att in names(P)){
  plotPerformanceOAT(P, sim_stoch, metric=att,plotType='base',attSel=attPerturb)
}

P.shuffle = calcPerformanceAttributes(clim=clim_ref,sim=sim_stoch_shuffle,attSel=attSel)

attSel = 'P_day_all_tot_dwellTime'
par(mfrow=c(5,5),mar=c(3,3,2,1),oma=c(1,1,3,1))
for (att in names(P.shuffle)){
  plotPerformanceOAT(P.shuffle, sim_stoch_shuffle, metric=att,plotType='base',attSel=attSel)
  title(attSel,outer=T)
}

####################

clim_sim = list(times=clim_ref$times,P=sim_stoch$Rep1$Target1$P$sim)
calculateAttributes(clim_sim,'P_day_all_P99')
calculateAttributes(clim_ref,'P_day_all_P99')

####################

attSel = c('P_day_all_tot','P_day_all_avgDSD','P_day_all_P99','P_day_all_nWet')
P = calcPerformanceAttributes(clim=clim_ref,sim=sim_stoch,attSel=attSel)
for (att in names(P)){
  plotPerformanceOAT(P, sim_stoch, metric=att,plotType='base',attSel=attPerturb)
}

####################

if (!is.null(Qobs)){
  
  source(paste0(runDirname,'GR4J_funcs.R'))
  
  Param = setup_cal_GR4J(dates,clim_ref$P,clim_ref$PET,Qobs)
  
  systemArgs = list(dates=dates,Param=Param)
  
  sysOutSim = runSystemModel(sim=sim_stoch,systemModel = GR4J_wrapper,systemArgs = systemArgs,metrics = c('meanQ','P99','P25'))
  
  sysOutClim = GR4J_wrapper(data = clim_ref, systemArgs = systemArgs,metrics = c('meanQ','P99','P25'))
  
  plotPerformanceOAT(performance = sysOutSim, sim=sim_stoch, metric = 'meanQ')
  plotPerformanceOAT(performance = sysOutSim, sim=sim_stoch, metric = 'P99')
  plotPerformanceOAT(performance = sysOutSim, sim=sim_stoch, metric = 'P25')
  
}

####################

dev.off()


