rm(list=ls())

devtools::load_all()

#library(foreSIGHT)

clim = foreSIGHT::convert_climYMD_POSIXct(tank_obs)

########################################

attsEval = c('P_day_all_tot_m','P_day_SON_tot_m','P_day_DJF_tot_m','P_day_MAM_tot_m','P_day_JJA_tot_m',
             'P_day_all_P99','P_day_SON_P99','P_day_DJF_P99','P_day_MAM_P99','P_day_JJA_P99')


# ########################################
# 
# modelSelection = list()
# modelSelection$modelType = list()
# modelSelection$modelType$P = "distScaling"
# modelSelection$modelParameterVariation = list()
# modelSelection$modelParameterVariation$P = "ann"
# modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
# controlFile = paste0(tempdir(), "\\eg_controlFile.json")
# write(modelSelectionJSON, file = controlFile)
# 
# ########################################
# 
# attPerturb = c("P_day_all_P99")
# attHold = c("P_day_all_tot_m")
# 
# attPerturbType = "regGrid"
# attPerturbSamp = c(7)
# attPerturbMin = c(0.8)
# attPerturbMax = c(1.4)
# 
# # create the exposure space
# expSpace = createExpSpace(attPerturb = attPerturb,
#                           attPerturbSamp = attPerturbSamp,
#                           attPerturbMin = attPerturbMin,
#                           attPerturbMax = attPerturbMax,
#                           attPerturbType = attPerturbType,
#                           attHold = attHold)
# 
# ########################################
# 
# sim = generateScenarios(reference = clim,
#                               expSpace = expSpace,
#                               controlFile = controlFile,
#                               seedID = 1)
# 
# plotScenarios(sim)
# 
# ########################################
# 
# # attSel = c('P_day_all_tot_m','P_day_SON_tot_m','P_day_DJF_tot_m','P_day_MAM_tot_m','P_day_JJA_tot_m',
# #              'P_day_all_P99','P_day_SON_P99','P_day_DJF_P99','P_day_MAM_P99','P_day_JJA_P99')
# # 
# # #attSel = unique(c(attPerturb,attHold,attOther))
# # # attSel = unique(c(attPerturb,attOther))
# # 
# # P = calcPerformanceAttributes(clim=convert_climYMD_POSIXct(tank_obs),sim=sim,attSel=attSel)
# # 
# # attSel='P_day_all_P99'
# # par(mfrow=c(2,5),mar=c(4,4,2,1))
# # for (att in names(P)){
# #   # plotPerformanceSpace(P, sim, metric=metric)
# #   plotPerformanceOAT(P, sim, metric=att,col='black',use_ggplot = F,attSel=attSel)
# # }
# 
# plot_attributes_vs_perturbed_OAT(clim,sim,attsEval,attPerturb='P_day_all_P99')

########################################







########################################

# modelSelection = list()
# modelSelection$modelType = list()
# modelSelection$modelType$P = "latent"
# modelSelection$modelParameterVariation = list()
# modelSelection$modelParameterVariation$P = "ann"
# modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
# controlFile = paste0(tempdir(), "\\eg_controlFile.json")
# write(modelSelectionJSON, file = controlFile)
# 
# ########################################
# 
# attPerturb = c("P_day_all_P99")
# attHold = c("P_day_all_tot_m","P_day_all_nWet_m","P_day_all_avgWSD_m")
# 
# attPerturbType = "regGrid"
# attPerturbSamp = c(7)
# attPerturbMin = c(0.8)
# attPerturbMax = c(1.4)
# 
# # create the exposure space
# expSpace = createExpSpace(attPerturb = attPerturb,
#                           attPerturbSamp = attPerturbSamp,
#                           attPerturbMin = attPerturbMin,
#                           attPerturbMax = attPerturbMax,
#                           attPerturbType = attPerturbType,
#                           attHold = attHold)
# 
# ########################################
# 
# sim = generateScenarios(reference = clim,
#                         expSpace = expSpace,
#                         controlFile = controlFile,
#                         seedID = 1)
# 
# plotScenarios(sim)

########################################

# attSel = c('P_day_all_tot_m','P_day_SON_tot_m','P_day_DJF_tot_m','P_day_MAM_tot_m','P_day_JJA_tot_m',
#            'P_day_all_P99','P_day_SON_P99','P_day_DJF_P99','P_day_MAM_P99','P_day_JJA_P99')
# 
# #attSel = unique(c(attPerturb,attHold,attOther))
# # attSel = unique(c(attPerturb,attOther))
# 
# P = calcPerformanceAttributes(clim=convert_climYMD_POSIXct(tank_obs),sim=sim,attSel=attSel)
# 
# attSel='P_day_all_P99'
# par(mfrow=c(2,5),mar=c(4,4,2,1))
# for (att in names(P)){
#   # plotPerformanceSpace(P, sim, metric=metric)
#   plotPerformanceOAT(P, sim, metric=att,col='black',use_ggplot = F,attSel=attSel)
# }

# plot_attributes_vs_perturbed_OAT(clim=clim,sim=sim,attsEval=attsEval,attPerturb='P_day_all_P99')








########################################

# modelSelection = list()
# modelSelection$modelType = list()
# modelSelection$modelType$P = "latent"
# modelSelection$modelParameterVariation = list()
# modelSelection$modelParameterVariation$P = "seas"
# modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
# controlFile = paste0(tempdir(), "\\eg_controlFile.json")
# write(modelSelectionJSON, file = controlFile)
# 
# 
# ########################################
# 
# attPerturb = c("P_day_all_P99")
# attHold = c("P_day_all_tot_m","P_day_all_nWet","P_day_all_avgWSD_m")
# 
# attPerturbType = "regGrid"
# attPerturbSamp = c(5)
# attPerturbMin = c(0.8)
# attPerturbMax = c(1.2)
# # attPerturbSamp = c(1)
# # attPerturbMin = c(0.8)
# # attPerturbMax = c(0.8)
# 
# # create the exposure space
# expSpace = createExpSpace(attPerturb = attPerturb,
#                           attPerturbSamp = attPerturbSamp,
#                           attPerturbMin = attPerturbMin,
#                           attPerturbMax = attPerturbMax,
#                           attPerturbType = attPerturbType,
#                           attHold = attHold)
# 
# attsAll = c(attPerturb,attHold)
# for (att in attsAll){
#   for (seas in c('DJF','MAM','JJA','SON')){
#     att.seas = gsub('all',seas,att)
#     expSpace$targetMat[att.seas] = expSpace$targetMat[att]
# #    expSpace$attPerturb = c(expSpace$attPerturb,att.seas)
#     expSpace$attHold = c(expSpace$attHold,att.seas)
#   }
# }
# #attsAll = c(expSpace$attPerturb,expSpace$attHold)
# #expSpace$targetMat = expSpace$targetMat[,attsAll]
# 
# sim = generateScenarios(reference = clim,
#                         expSpace = expSpace,
#                         controlFile = controlFile,
#                         seedID = 1,numReplicates = 5)
# 
# plotScenarios(sim)
# 
# ########################################
# 
# plot_attributes_vs_perturbed_OAT(clim,sim,attsEval,attPerturb='P_day_all_P99')
  
  
# P = calcPerformanceAttributes(clim=clim,sim=sim,attSel=attSel)
# 
# attSel='P_day_all_P99'
# par(mfrow=c(2,5),mar=c(4,4,2,1))
# for (att in names(P)){
#   plotPerformanceOAT(P, sim, metric=att,plotType='base',attSel=attSel)
# }











########################################

modelSelection = list()
modelSelection$modelType = list()
modelSelection$modelType$P = "latent"
modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "har"
modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")
write(modelSelectionJSON, file = controlFile)

########################################

attPerturb = c("P_day_all_P99")
attHold = c("P_day_all_tot_m","P_day_all_nWet","P_day_all_avgWSD_m")

attPerturbType = "regGrid"
attPerturbSamp = c(3)
attPerturbMin = c(0.8)
attPerturbMax = c(1.2)
# attPerturbSamp = c(1)
# attPerturbMin = c(0.8)
# attPerturbMax = c(0.8)

# create the exposure space
expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType,
                          attHold = attHold,
                          attTied = c(attPerturb,attHold),
                          tieType = 'seas')

# attsAll = c(attPerturb,attHold)
# for (att in attsAll){
#   for (seas in c('DJF','MAM','JJA','SON')){
#     att.seas = gsub('all',seas,att)
#     expSpace$targetMat[att.seas] = expSpace$targetMat[att]
#     expSpace$attTied = c(expSpace$attHold,att.seas)
#   }
# }

time1 = Sys.time()
simHar = generateScenarios(reference = clim,
                           expSpace = expSpace,
                           controlFile = controlFile,
                           seedID = 1,numReplicates = 2)
time2 = Sys.time()
print(time2-time1)

# time1 = Sys.time()
# simHar = generateScenarios(reference = clim,
#                         expSpace = expSpace,
#                         controlFile = controlFile,
#                         seedID = 1,numReplicates = 2,cores=5)
# time2 = Sys.time()
# print(time2-time1)
# 
# time1a = Sys.time()
# simHar = generateScenarios(reference = clim,
#                            expSpace = expSpace,
#                            controlFile = controlFile,
#                            seedID = 1,numReplicates = 2,cores=1)
# time2a = Sys.time()
# print(time2a-time1a)

plotScenarios(simHar)

########################################

attsEval = c('P_day_all_tot_m','P_day_SON_tot_m','P_day_DJF_tot_m','P_day_MAM_tot_m','P_day_JJA_tot_m',
             'P_day_all_P99','P_day_SON_P99','P_day_DJF_P99','P_day_MAM_P99','P_day_JJA_P99',
             'P_day_all_nWet_m','P_day_SON_nWet_m','P_day_DJF_nWet_m','P_day_MAM_nWet_m','P_day_JJA_nWet_m',
             'P_day_all_avgWSD','P_day_SON_avgWSD','P_day_DJF_avgWSD','P_day_MAM_avgWSD','P_day_JJA_avgWSD')

par(mfrow=c(4,5))
plot_attributes_vs_perturbed_OAT(clim,simHar,attsEval,attPerturb='P_day_all_P99',
                                 baseSettings=list(bias_base_thresh = 5,slope_thresh = 1.5))


