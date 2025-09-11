# script for generating single-site perturbed rainfall and PET data for Loxton for Source w crop model 

# clear all variables
rm(list=ls())

#----------------------------------------------------------------
# Set paths and get data
#----------------------------------------------------------------

# path to foreSIGHT code from repo  
foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT/'
#foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT_current_github/foreSIGHT/'

# path for data files
dataDir = 'C:/Users/a1065639/Box/2025_DEW_foreSIGHT/Data/'

# load foreSIGHT package
setwd(foreSIGHTDir)
devtools::load_all()

# remove.packages("foreSIGHT")
# .rs.restartR()
# devtools::install()
#library(foreSIGHT)

# climate reference data file name
fname = paste0(dataDir,'24007_LOXTON_Dec2023.csv')

# read in evap reference data
d = read.csv(fname,header=T,skip=4017) # note there is different format for dates pre and post 1900. have chosen to not read in pre-1900 data.

date = as.Date(d[,1],format = '%d/%m/%Y')
P = d[,2]
PET = d[,3]

# select time periods of data to use for simulations
startYr = 1976
endYr = 2005
# startYr = 1980
# endYr = 1989
year = as.integer(format(date,'%Y'))
keep = which((year>=startYr)&(year<=endYr))
# create reference data 
clim_ref <- list(times = as.POSIXct(date[keep],tz='UTC'),
                 P = P[keep],
                 PET = PET[keep])

#-----------------------------------------------------------------
# Step A: Create an exposure space
#-----------------------------------------------------------------

# # specify perturbed attributes
# attPerturb <- c("PET_day_all_avg_m", "P_day_all_tot_m") 
# 
# # specify perturbation type and minimum-maximum ranges of the perturbed attributes
# attPerturbType <- "regGrid"
# attPerturbSamp <- c(7,7)     
# attPerturbMin <- c(1, 0.7)
# attPerturbMax <- c(1.3, 1.3)
# # create the exposure space using foreSIGHT
# expSpace <- createExpSpace(attPerturb = attPerturb,
#                            attPerturbSamp = attPerturbSamp,            # set to null
#                            attPerturbMin = attPerturbMin,
#                            attPerturbMax = attPerturbMax,
#                            attPerturbType = attPerturbType,
#                            attHold = NULL)
# 

#-----------------------------------------------------------------
# Step B: Generate perturbed time series
#-----------------------------------------------------------------



# generate simulations using simple scaling with seasonality 
# sim_scaling <- generateScenarios(reference = clim_ref,             # input observed data
#                          expSpace = expSpace,        # exposure space created by the user
#                          controlFile = "scaling")    # using simple scaling with seasonality



######################################################################

# mvFunc_cor = function(data.1,data.2){
#   return(cor(data.1,data.2,use='pairwise.complete.obs'))
# }
# 
# mvFunc_avgWetDay = function(data.1,data.2){
#   return(mean(data.1[data.2>0],na.rm=T))
# }
# 
# mvFunc_sdWetDay = function(data.1,data.2){
#   return(sd(data.1[data.2>0],na.rm=T))
# }
# 
# mvFunc_avgDryDay = function(data.1,data.2){
#   return(mean(data.1[data.2==0],na.rm=T))
# }
# 
# mvFunc_sdDryDay = function(data.1,data.2){
#   return(sd(data.1[data.2==0],na.rm=T))
# }
# 
# func_sd = function(data){
#   return(sd(data,na.rm=T))
# }
# 
# func_xP90 = function(data){
#   P90 = quantile(data,probs = 0.9,na.rm=T,names=F)
#   return(P90)
# }
# 
# mvFunc_xP90WetDay = function(data.1,data.2){
#   return(quantile(data.1[data.2>0],probs=0.9,na.rm=T,names=F))
# }
# 
# mvFunc_xP90DryDay = function(data.1,data.2){
#   return(quantile(data.1[data.2==0],probs=0.9,na.rm=T,names=F))
# }
# 
# func_cv = function(data){
#   m = mean(data,na.rm=T)
#   if (m==0){
#     cv = 9999.
#   } else {
#     cv = sd(data,na.rm=T)/m
#   }
#   return(cv)
# }
# 
# mvFunc_cvWetDay = function(data.1,data.2){
#   return(func_cv(data.1[data.2>0]))
# }
# 
# mvFunc_cvDryDay = function(data.1,data.2){
#   return(func_cv(data.1[data.2==0]))
# }



modelSelection = list()

modelSelection$modelType = list()
modelSelection$modelType$P = "latent"
modelSelection$modelType$PET = "wgenLM"

modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "seas"
# modelSelection$modelParameterVariation$P = "har"
# modelSelection$modelParameterVariation$PET = "seasWD"
modelSelection$modelParameterVariation$PET = "harWD"

# modelSelection[["modelParameterBounds"]] <- list()
# modelSelection[["modelParameterBounds"]][["P"]] <- list()
# modelSelection[["modelParameterBounds"]][["P"]][["sigma.SON"]] <- c(0.001, 20)
# modelSelection[["modelParameterBounds"]][["P"]][["sigma.DJF"]] <- c(0.001, 20)
# modelSelection[["modelParameterBounds"]][["P"]][["sigma.MAM"]] <- c(0.001, 20)
# modelSelection[["modelParameterBounds"]][["P"]][["sigma.JJA"]] <- c(0.001, 20)
# modelSelection[["modelParameterBounds"]][["P"]][["mu.SON"]] <- c(-20, 1)
# modelSelection[["modelParameterBounds"]][["P"]][["mu.DJF"]] <- c(-20, 1)
# modelSelection[["modelParameterBounds"]][["P"]][["mu.MAM"]] <- c(-20, 1)
# modelSelection[["modelParameterBounds"]][["P"]][["mu.JJA"]] <- c(-20, 1)
# modelSelection[["modelParameterBounds"]][["P"]][["lambda.SON"]] <- c(0.5, 5)
# modelSelection[["modelParameterBounds"]][["P"]][["lambda.DJF"]] <- c(0.5, 5)
# modelSelection[["modelParameterBounds"]][["P"]][["lambda.MAM"]] <- c(0.5, 5)
# modelSelection[["modelParameterBounds"]][["P"]][["lambda.JJA"]] <- c(0.5, 5)

modelSelection[["modelParameterBounds"]] <- list()
modelSelection[["modelParameterBounds"]][["P"]] <- list()
modelSelection[["modelParameterBounds"]][["P"]][["alpha.SON"]] <- c(0.1, 0.9)
modelSelection[["modelParameterBounds"]][["P"]][["alpha.DJF"]] <- c(0.1, 0.9)
modelSelection[["modelParameterBounds"]][["P"]][["alpha.MAM"]] <- c(0.1, 0.9)
modelSelection[["modelParameterBounds"]][["P"]][["alpha.JJA"]] <- c(0.1, 0.9)
modelSelection[["modelParameterBounds"]][["P"]][["sigma.SON"]] <- c(0.5, 8)
modelSelection[["modelParameterBounds"]][["P"]][["sigma.DJF"]] <- c(0.5, 8)
modelSelection[["modelParameterBounds"]][["P"]][["sigma.MAM"]] <- c(0.5, 8)
modelSelection[["modelParameterBounds"]][["P"]][["sigma.JJA"]] <- c(0.5, 8)
modelSelection[["modelParameterBounds"]][["P"]][["mu.SON"]] <- c(-8, 1)
modelSelection[["modelParameterBounds"]][["P"]][["mu.DJF"]] <- c(-8, 1)
modelSelection[["modelParameterBounds"]][["P"]][["mu.MAM"]] <- c(-8, 1)
modelSelection[["modelParameterBounds"]][["P"]][["mu.JJA"]] <- c(-8, 1)
modelSelection[["modelParameterBounds"]][["P"]][["lambda.SON"]] <- c(1, 5)
modelSelection[["modelParameterBounds"]][["P"]][["lambda.DJF"]] <- c(1, 5)
modelSelection[["modelParameterBounds"]][["P"]][["lambda.MAM"]] <- c(1, 5)
modelSelection[["modelParameterBounds"]][["P"]][["lambda.JJA"]] <- c(1, 5)

# 
# modelSelection[["modelParameterBounds"]][["PET"]] <- list()
# modelSelection[["modelParameterBounds"]][["PET"]][["cor0.SON"]] <- c(0.1, 0.99)
# modelSelection[["modelParameterBounds"]][["PET"]][["cor0.DJF"]] <- c(0.1, 0.99)
# modelSelection[["modelParameterBounds"]][["PET"]][["cor0.MAM"]] <- c(0.1, 0.99)
# modelSelection[["modelParameterBounds"]][["PET"]][["cor0.JJA"]] <- c(0.1, 0.99)
# modelSelection[["modelParameterBounds"]][["PET"]][["muW.SON"]] <- c(0, 10)
# modelSelection[["modelParameterBounds"]][["PET"]][["muW.DJF"]] <- c(0, 10)
# modelSelection[["modelParameterBounds"]][["PET"]][["muW.MAM"]] <- c(0, 10)
# modelSelection[["modelParameterBounds"]][["PET"]][["muW.JJA"]] <- c(0, 10)
# modelSelection[["modelParameterBounds"]][["PET"]][["sigmaW.SON"]] <- c(0.01, 5)
# modelSelection[["modelParameterBounds"]][["PET"]][["sigmaW.DJF"]] <- c(0.01, 5)
# modelSelection[["modelParameterBounds"]][["PET"]][["sigmaW.MAM"]] <- c(0.01, 5)
# modelSelection[["modelParameterBounds"]][["PET"]][["sigmaW.JJA"]] <- c(0.01, 5)
# modelSelection[["modelParameterBounds"]][["PET"]][["muD.SON"]] <- c(0, 15)
# modelSelection[["modelParameterBounds"]][["PET"]][["muD.DJF"]] <- c(0, 15)
# modelSelection[["modelParameterBounds"]][["PET"]][["muD.MAM"]] <- c(0, 15)
# modelSelection[["modelParameterBounds"]][["PET"]][["muD.JJA"]] <- c(0, 15)
# modelSelection[["modelParameterBounds"]][["PET"]][["sigmaD.SON"]] <- c(0.01, 5)
# modelSelection[["modelParameterBounds"]][["PET"]][["sigmaD.DJF"]] <- c(0.01, 5)
# modelSelection[["modelParameterBounds"]][["PET"]][["sigmaD.MAM"]] <- c(0.01, 5)
# modelSelection[["modelParameterBounds"]][["PET"]][["sigmaD.JJA"]] <- c(0.01, 5)

# parNam=c("alpha.SON","alpha.DJF","alpha.MAM","alpha.JJA",
#          "sigma.SON","sigma.DJF","sigma.MAM","sigma.JJA",
#          "mu.SON","mu.DJF","mu.MAM","mu.JJA",
#          "lambda.SON","lambda.DJF","lambda.MAM","lambda.JJA"),
# 
# minBound=c(0,0,0,0,
#            0.001,0.001,0.001,0.001,
#            -15,-15,-15,-15,
#            1,1,1,1),
# maxBound=c(0.999,0.999,0.999,0.999,
#            10,10,10,10,
#            1,1,1,1,
#            4,4,4,4)

modelSelection[["optimisationArguments"]] <- list()
modelSelection[["optimisationArguments"]][["nMultiStart"]] <- 5
modelSelection[["optimisationArguments"]][["OFtol"]] <- 0.02
# modelSelection[["optimisationArguments"]][['RGN.control']]=list(iterMax=50)

#modelSelection[["penaltyAttributes"]] <- c("P_day_all_tot", "P_day_all_xP99overPave",
#                                           "P_day_all_nWet", "P_day_all_avgDSD")#,

#modelSelection[["penaltyAttributes"]] <- c("P_day_all_tot", "P_day_all_xP99overPave",
#                                             "P_day_all_nWet", "P_day_all_avgDSD",
#                                           "PET_day_all_avg","PET_day_all_cor","PET_day_all_cv")#,"P_year_all_cv")
modelSelection[["penaltyAttributes"]] <- c("P_day_all_tot", "P_day_all_P99",
                                           "P_day_all_nWet", "P_day_all_avgDSD",
                                           "PET_day_all_avg","PET_day_all_cor","PET_day_all_cv")#,"P_year_all_cv")
#modelSelection[["penaltyWeights"]] <-c(2,2,2,2)
#modelSelection[["penaltyWeights"]] = rep(3,8)
modelSelection[["penaltyWeights"]] = rep(3,7)
#modelSelection[["penaltyWeights"]] = rep(3,4)
# modelSelection[["penaltyWeights"]] = rep(0,4)

# modelSelection$postProcessing = list()
# modelSelection$postProcessing$P = list()
# modelSelection$postProcessing$P$types = c('annVar','annCor')

modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")
write(modelSelectionJSON, file = controlFile)

####################

# #attPerturb = c("P_day_all_tot_m")
# #attHold = c('P_day_all_P99','P_day_all_nWet','P_day_all_avgDSD_m')
# 
# attPerturb = c('P_day_all_avgDSD_m')
# #attHold = c('P_day_all_tot_m','P_day_all_P99','P_day_all_nWet',
# #            "PET_day_all_cor","mv.PET.P_day_all_meanWetDay","mv.PET.P_day_all_sdWetDay",
# #            "mv.PET.P_day_all_meanDryDay","mv.PET.P_day_all_sdDryDay")
# attHold = c('P_day_all_tot_m','P_day_all_P99','P_day_all_nWet',
#             "PET_day_all_cor","PET_day_all_mean","PETday_all_sd")
# 
# attPerturbType = "regGrid"
# attPerturbSamp = c(2)
# attPerturbMin = c(1)
# attPerturbMax = c(1.5)


# attHold = c('P_day_all_P99','P_day_all_nWet','P_day_all_avgDSD_m',
#             "P_day_SON_tot_m",'P_day_SON_P99','P_day_SON_nWet','P_day_SON_avgDSD_m',
#             'P_day_DJF_tot_m','P_day_DJF_P99','P_day_DJF_nWet','P_day_DJF_avgDSD_m',
#             'P_day_MAM_tot_m','P_day_MAM_P99','P_day_MAM_nWet','P_day_MAM_avgDSD_m',
#             'P_day_JJA_tot_m','P_day_JJA_P99','P_day_JJA_nWet','P_day_JJA_avgDSD_m')#,
#             # "PET_day_all_cor","mv.PET.P_day_all_meanWetDay","mv.PET.P_day_all_sdWetDay",
#             # "mv.PET.P_day_all_meanDryDay","mv.PET.P_day_all_sdDryDay",
#             # "PET_day_SON_cor","mv.PET.P_day_SON_meanWetDay","mv.PET.P_day_SON_sdWetDay",
#             # "mv.PET.P_day_SON_meanDryDay","mv.PET.P_day_SON_sdDryDay",
#             # "PET_day_DJF_cor","mv.PET.P_day_DJF_meanWetDay","mv.PET.P_day_DJF_sdWetDay",
#             # "mv.PET.P_day_DJF_meanDryDay","mv.PET.P_day_DJF_sdDryDay",
#             # "PET_day_MAM_cor","mv.PET.P_day_MAM_meanWetDay","mv.PET.P_day_MAM_sdWetDay",
#             # "mv.PET.P_day_MAM_meanDryDay","mv.PET.P_day_MAM_sdDryDay",
#             # "PET_day_JJA_cor","mv.PET.P_day_JJA_meanWetDay","mv.PET.P_day_JJA_sdWetDay",
#             # "mv.PET.P_day_JJA_meanDryDay","mv.PET.P_day_JJA_sdDryDay"
#             # )


# attPerturb = c('P_day_all_avgDSD_m',"PET_day_all_avg")
# attHold = c('P_day_all_tot_m','P_day_all_P99','P_day_all_nWet',
#             "PET_day_all_cor","PET_day_all_xP90")

#attPerturbType = "regGrid"
#attPerturbSamp = c(2,2)
#attPerturbMin = c(1,1)
#attPerturbMax = c(1.5,1.5)

# attPerturbType = "regGrid"
# attPerturbSamp = c(1,1)
# attPerturbMin = c(1.3,1.3)
# attPerturbMax = c(1.3,1.3)

#attPerturb = c('P_day_all_avgDSD_m')
#attHold = c('P_day_all_tot_m','P_day_all_P99','P_day_all_nWet',
#            "PET_day_all_avg","PET_day_all_cor","PET_day_all_xP90")

#########################

attPerturb = c('P_day_all_tot')
attHold = c('P_day_all_avgDSD','P_day_all_P99','P_day_all_nWet',
            "PET_day_all_avg","PET_day_all_cor","PET_day_all_xP90")

attPerturbType = "regGrid"
attPerturbSamp = c(5)
attPerturbMin = c(0.7)
attPerturbMax = c(1.1)

#########################

# attPerturb = c('P_day_all_avgDSD')
# attHold = c('P_day_all_tot','P_day_all_xP99overPave','P_day_all_nWet',
#             "PET_day_all_avg","PET_day_all_cor","PET_day_all_cv")
# 
# attPerturbType = "regGrid"
# attPerturbSamp = c(5)
# attPerturbMin = c(1.)
# attPerturbMax = c(1.4)

##########################

# attPerturb = c("PET_day_all_avg")
# 
# #attHold = c('P_day_all_tot','P_day_all_avgDSD','P_day_all_xP99overPave','P_day_all_nWet',
# #            "PET_day_all_cor","PET_day_all_cv")
# 
# attHold = c('P_day_all_tot','P_day_all_avgDSD','P_day_all_P99','P_day_all_nWet',
#             "PET_day_all_cor","PET_day_all_cv")
# 
# attPerturbType = "regGrid"
# #attPerturbSamp = c(5)
# #attPerturbMin = c(1)
# #attPerturbMax = c(1.2)
# 
# attPerturbSamp = c(1)
# attPerturbMin = c(1)
# attPerturbMax = c(1)

#########################

# attPerturb = c('P_day_all_tot','P_year_all_cor')
# 
# attHold = c("PET_day_all_avg",'P_day_all_avgDSD','P_day_all_xP99overPave','P_day_all_nWet',
#             "PET_day_all_cor","PET_day_all_cv",'P_year_all_cv')
# 
# #attPerturbType = "regGrid"
# #attPerturbSamp = c(5)
# #attPerturbMin = c(0.8)
# #attPerturbMax = c(1.2)
# 
# attPerturbType = "regGrid"
# attPerturbSamp = c(1,1)
# attPerturbMin = c(1,0)
# attPerturbMax = c(1,0)

#########################

# attPerturbType = "regGrid"
# attPerturbSamp = c(5)
# attPerturbMin = c(1)
# attPerturbMax = c(1.4)
# 
# # attPerturbSamp = c(1)
# # attPerturbMin = c(0.7)
# # attPerturbMax = c(0.7)

#########################

# attPerturbSamp = c(1)
# attPerturbMin = c(1)
# attPerturbMax = c(1)

#########################

#attPerturb = c('P_day_all_tot')
#attHold = c('P_day_all_avgDSD','P_day_all_P99','P_day_all_nWet')

# attPerturb = c('P_day_all_tot')
# attHold = c('P_day_all_avgDSD','P_day_all_xP99overPave','P_day_all_nWet')

# attPerturbType = "regGrid"
# attPerturbSamp = c(5)
# attPerturbMin = c(0.7)
# attPerturbMax = c(1.1)

# attPerturbType = "regGrid"
# attPerturbSamp = c(1)
# attPerturbMin = c(0.8)
# attPerturbMax = c(0.8)


# attPerturbType = "regGrid"
# attPerturbSamp = c(5)
# attPerturbMin = c(0.8)
# attPerturbMax = c(1.2)


#########################

#attPerturb = c('P_day_all_avgDSD')
#attHold = c('P_day_all_tot','P_day_all_xP99overPave','P_day_all_nWet')

# attPerturb = c('P_day_all_avgDSD')
# attHold = c('P_day_all_tot','P_day_all_P99','P_day_all_nWet')
# 
# # attPerturbType = "regGrid"
# # attPerturbSamp = c(5)
# # attPerturbMin = c(0.7)
# # attPerturbMax = c(1.1)
# 
# attPerturbType = "regGrid"
# attPerturbSamp = c(1)
# attPerturbMin = c(1.3)
# attPerturbMax = c(1.3)

#########################

# create the exposure space
expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType,
                          attHold = attHold,
                          targetTypes = list(P_year_all_cor='val') 
                          )#,
#                          attTied = c(attPerturb,attHold),
#                          tieType = 'seas')

pause

# attTied = list(wDdD=c("PET_day_all_avg","PET_day_all_cv"),
#               seas='allTargets')
# attTied = list(seas='allTargets')

attTied = list(wDdD=c("PET_day_all_avg","PET_day_all_cv"))
expSpace = tieAttributes(expSpace,attTied)

atts = colnames(expSpace$targetMat)
seasTied = atts[!atts%in%c('P_year_all_cv','P_year_all_cor')]
attTied = list(seas=seasTied)
expSpace = tieAttributes(expSpace,attTied)

i=which(grepl(pattern = 'xP99overPave',x = colnames(expSpace$targetMat)))

expSpace$targetMat[,i] = 1/expSpace$targetMat[,'P_day_all_tot']

time.1 = Sys.time()
sim_stoch = generateScenarios(reference = clim_ref,
                              expSpace = expSpace,
                              controlFile = controlFile,
                              seedID = 1,
                              numReplicates = 1,
                              cores = 1)
time.2 = Sys.time()
print(time.2-time.1)

fname = paste0('testing_2.0/Loxton_',attPerturb,'.RData') 
save(file=fname,sim_stoch,clim_ref)

plotScenarios(sim_stoch)

attSel = colnames(expSpace$targetMat)

P = calcPerformanceAttributes(clim=clim_ref,sim=sim_stoch,attSel=attSel)

par(mfrow=c(5,5),mar=c(4,5,2,1))
for (att in names(P)){
  plotPerformanceOAT(P, sim_stoch, metric=att,plotType='base',attSel=attPerturb)
}

####################

clim_sim = list(times=clim_ref$times,P=sim_stoch$Rep5$Target1$P$sim)
calculateAttributes(clim_sim,'P_day_all_P99')
calculateAttributes(clim_ref,'P_day_all_P99')

####################

attSel = c('P_day_all_tot','P_day_all_avgDSD','P_day_all_P99','P_day_all_nWet')
P = calcPerformanceAttributes(clim=clim_ref,sim=sim_stoch,attSel=attSel)
for (att in names(P)){
  plotPerformanceOAT(P, sim_stoch, metric=att,plotType='base',attSel=attPerturb)
}


####################

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
