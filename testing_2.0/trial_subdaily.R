rm(list=ls())

devtools::load_all()

#devtools::load_all('C:/Users/a1065639/Work/RGN/')

# load('testing_2.0/AWS_023034_accumHandling.spread_processed_aggPeriod.1hour_aggDataThresh.0.8_agg.RData')
# 
# agg_data$date
# agg_data$P
# 
# clim = list(year=as.integer(format(agg_data$date,'%Y')),
#             month=as.integer(format(agg_data$date,'%m')),
#             day=as.integer(format(agg_data$date,'%d')),
#             hour=as.integer(format(agg_data$date,'%H')),
#             P=agg_data$P,
#             timeStep='1 hour')
# 
# calculateAttributes(clim,'P_ann_tot_m')


library('BLRPM')

timeStart = as.POSIXct('2000/01/01 00:00:00',tz='UTC')
timeEnd = as.POSIXct('2000/12/31 23:00:00',tz='UTC')

# timeStart = as.POSIXct('2000/01/01 00:00:00',tz='UTC')
# timeEnd = as.POSIXct('2009/12/31 23:00:00',tz='UTC')

times = seq(timeStart,timeEnd,by='hours')

nTimes = length(times)

lambda = 0.02
gamma = 1/10
beta = 0.3
eta = 2
mux = 4
t.sim = nTimes
# t.acc = t.sim
# interval = 1
# offset = 0
seed = 1 
set.seed(seed)
# simulation0 = BLRPM.sim(lambda,gamma,beta,eta,mux,t.sim)
# 
# BLRPM::BLRPM.est(simulation0$RR)

simulation = SWGsim.BLRPM(SWGpar=list(lambda=lambda,gamma=gamma,beta=beta,eta=eta,mux=mux),
                           nTimes=nTimes,
                           randomTerm=list(seed=seed))

clim = list(times=times,
            P=simulation)
########################################

# func_var = function(data){
#   var(data)
# }
# func_cov = function(data){
#   a = acf(data,lag.max=1,type="covariance",plot=F)$acf[2,1,1]
#   return(a)
# } 
# func_probZero = function(data) length(data[data==0])/length(data)  

########################################

modelSelection = list()
modelSelection$modelType = list()
modelSelection$modelType$P = "BLRPM"
modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "ann"

# modelSelection[["modelParameterBounds"]] <- list()
# modelSelection[["modelParameterBounds"]][["P"]] <- list()
# modelSelection[["modelParameterBounds"]][["P"]][["lambda"]] <- c(0.01, 0.03)
# modelSelection[["modelParameterBounds"]][["P"]][["gamma"]] <- c(0.05, 0.2)
# modelSelection[["modelParameterBounds"]][["P"]][["beta"]] <- c(0.2, 0.4)
# modelSelection[["modelParameterBounds"]][["P"]][["eta"]] <- c(1.5, 2.5)
# modelSelection[["modelParameterBounds"]][["P"]][["mux"]] <- c(3, 5)

# modelSelection[["penaltyAttributes"]] <- c('P_day_all_tot')
# modelSelection[["penaltyWeights"]] <- c(10)

#modelSelection[["optimisationArguments"]] <- list()
#modelSelection[["optimisationArguments"]]$suggestions = c(0.02,0.1,0.3,2.0,4)
#modelSelection[["optimisationArguments"]][["optimizer"]] <- 'NM'
# modelSelection[["optimisationArguments"]][["optimizer"]] <- 'SCE'
#modelSelection[["optimisationArguments"]][["nMultiStart"]] <- 5

modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)

controlFile = paste0(tempdir(), "\\eg_controlFile.json")
write(modelSelectionJSON, file = controlFile)


########################################

# attPerturb = c('P_day_all_tot')
# attHold = c('P_hour_all_var','P_hour_all_cov','P_hour_all_probZero',
#            'P_3hour_all_var','P_3hour_all_cov','P_3hour_all_probZero',
#            'P_12hour_all_var','P_12hour_all_cov','P_12hour_all_probZero',
#            'P_day_all_var','P_day_all_cov','P_day_all_probZero')

attPerturb = c('P_day_all_tot')
attHold = c('P_hour_all_sd','P_hour_all_cor','P_hour_all_nWet',
            'P_3hour_all_sd','P_3hour_all_cor','P_3hour_all_nWet',
            'P_12hour_all_sd','P_12hour_all_cor','P_12hour_all_nWet',
            'P_day_all_sd','P_day_all_cor','P_day_all_nWet')


#attHold = c('P_hour_all_var','P_hour_all_cov','P_hour_all_probZero',
#            'P_3hour_all_var','P_3hour_all_cov','P_3hour_all_probZero',
#            'P_day_all_var','P_day_all_cov','P_day_all_probZero')

# attHold = c('P_hour_all_var','P_hour_all_cov','P_hour_all_probZero',
#            'P_day_all_var','P_day_all_cov','P_day_all_probZero')

# attHold = c('P_hour_all_var')

attPerturbType = "regGrid"
attPerturbSamp = c(1)
attPerturbMin = c(1.)
attPerturbMax = c(1.)

# create the exposure space
expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType,
                          attHold = attHold)

########################################

sim = generateScenarios(reference = clim,
                        expSpace = expSpace,
                        controlFile = controlFile,
                        seedID = 1)

plotScenarios(sim)


