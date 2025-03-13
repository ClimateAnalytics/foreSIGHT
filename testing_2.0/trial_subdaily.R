rm(list=ls())

devtools::load_all()

devtools::load_all('C:/Users/a1065639/Work/RGN/')

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

# timeStart = as.POSIXct('2000/01/01 00:00:00')
# timeEnd = as.POSIXct('2000/12/31 23:00:00')

timeStart = as.POSIXct('2000/01/01 00:00:00')
timeEnd = as.POSIXct('2009/12/31 23:00:00')

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

clim = list(year=as.integer(format(times,'%Y')),
            month=as.integer(format(times,'%m')),
            day=as.integer(format(times,'%d')),
            hour=as.integer(format(times,'%H')),
            P=simulation)

# calculateAttributes(clim,'P_ann_tot_m')
# 
# calculateAttributes(clim,'P_ann_P99.9')
# 
# calculateAttributes(clim,'P_ann_avgDSD')
# 
# calculateAttributes(clim,'P_ann_nWet')

########################################

func_max = function(data) max(data)

#atts = c('P_hour_all_max','P_day_all_P99','P_month_all_P99')
# atts = c('P_3hour_all_max','P_day_all_P99','P_month_all_P99')
# a=calculateAttributes(clim,atts)

########################################

func_var = function(data){
  # browser()
  var(data)
}
func_cov = function(data){
  if (is.na(sum(data))){browser()}
  a = acf(data,lag.max=1,type="covariance",plot=F)$acf[2,1,1]
  #print(a)
  return(a)
} 
func_probZero = function(data) length(data[data==0])/length(data)  

# atts = c('P_hour_all_tot',
#          'P_hour_all_var','P_hour_all_cov','P_hour_all_probZero',
#          'P_3hour_all_var','P_3hour_all_cov','P_3hour_all_probZero',
#          'P_12hour_all_var','P_12hour_all_cov','P_12hour_all_probZero',
#          'P_day_all_var','P_day_all_cov','P_day_all_probZero')
# 
# a=calculateAttributes(clim,atts)

########################################

modelSelection = list()
modelSelection$modelType = list()
modelSelection$modelType$P = "BLRPM"
modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "ann"

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

attPerturb = c('P_day_all_tot')
attHold = c('P_hour_all_var','P_hour_all_cov','P_hour_all_probZero',
           'P_3hour_all_var','P_3hour_all_cov','P_3hour_all_probZero',
           'P_12hour_all_var','P_12hour_all_cov','P_12hour_all_probZero',
           'P_day_all_var','P_day_all_cov','P_day_all_probZero')

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


