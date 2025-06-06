# #################################

#' @include default_parameters.R

modelInfoList[['P-ann-wgen']] = list(simVar='P',
                                        timeStep = '1 day',
                                        simPriority=1,
                                        npars=4,
                                        parNam = c('pdd','pwd','alpha','beta'),
                                     minBound=c(0.427, 0.088, 0.313, 0.043), 
                                     maxBound=c(0.95, 0.824, 0.998, 25.46))

                                       # minBound=c(0.427, 0.088, 0.313, 0.043), 
                                       # maxBound=c(0.998, 0.824, 0.998, 25.46))

# modelInfoList[['P-annDelta-wgen']] = list(simVar='P',
#                                             timeStep = '1 day',
#                                             simPriority=1,
#                                             npars=4,
#                                             parNam = c('delta.pdd','delta.pwd','delta.alpha','delta.beta'),
#                                             minBound=c(-1, -1, -1, -10),
#                                             maxBound=c(1, 1, 1, 10.))

modelInfoList[['P-seas-wgen']] = list(simVar='P',
                                        timeStep = '1 day',
                                        simPriority=1,
                                        npars=16,
                                        parNam = c('pdd.SON','pdd.DJF','pdd.MAM','pdd.JJA',
                                                     'pwd.SON','pwd.DJF','pwd.MAM','pwd.JJA',
                                                     'alpha.SON','alpha.DJF','alpha.MAM','alpha.JJA',
                                                     'beta.SON','beta.DJF','beta.MAM','beta.JJA'),
                                        minBound=c(0.389, 0.334, 0.375, 0.277,
                                                   0.078, 0.079, 0.084, 0.036,
                                                   0.295,	0.303, 0.309,	0.257,
                                                   0.043,	0.046, 0.048, 0.034), #Aus 3stdev hard bounds
                                        maxBound=c(0.997, 0.989, 0.994, 0.998,
                                                   0.85, 0.714, 0.714, 0.808,
                                                   0.998, 0.998, 0.998, 0.998,
                                                   15.716, 30.08, 27.877, 21.193))

###########################

modelInfoList[['P-har-wgen']] = list(simVar='P',
                                        timeStep = '1 day',
                                        simPriority=1,
                                        npars=12,
                                        parNam = c('pdd.m','pdd.amp','pdd.ang',
                                                     'pwd.m','pwd.amp','pwd.ang',
                                                     'alpha.m','alpha.amp','alpha.ang',
                                                     'beta.m','beta.amp','beta.ang'),
                                        minBound=c(0.476, 0.006, 0,
                                                   0.093, 0.004, 0,
                                                   0.33, 0.002, 0,
                                                   0.085, 0.028, 0),
                                        maxBound=c(0.950, 0.557, 6.28,
                                                   0.728, 0.519, 6.28,
                                                   0.950, 0.600, 6.28,
                                                   15.00, 10, 6.28))



#################################
 
parManager.wgen = function(parS, SWGparameterization, datInd, auxInfo=NULL){
  
  parTS = list()
  
  nTimes = datInd$nTimes
  
  # basePar = c(0.8004759,0.3494826,0.4233287,7.6178369)
  # names(basePar) = c('pdd','pwd','alpha','beta')
  if (SWGparameterization=='ann'){
    for (par in names(parTS)){
      parTS[[par]] = rep(parS[par],nTimes)
    }
#    parTS$pdd <- rep(parS['pdd'],datInd$nTimes)
#    parTS$pwd <- rep(parS['pwd'],datInd$nTimes)
#    parTS$alpha <- rep(parS['alpha'],datInd$nTimes)
#    parTS$beta <- rep(parS['beta'],datInd$nTimes)
  # } else if (SWGparameterization=='annDelta'){
  #   pdd = basePar['pdd']+parS['delta.pdd']#; pdd = max(min(pdd,1),0)
  #   pwd = basePar['pwd']+parS['delta.pwd']#; pwd = max(min(pwd,1),0)
  #   alpha = basePar['alpha']+parS['delta.alpha']#; alpha = max(alpha,0)
  #   beta = basePar['beta']+parS['delta.beta']#; beta = max(beta,0)
  #   parTS$pdd <- rep(pdd,nTimes)
  #   parTS$pwd <- rep(pwd,nTimes)
  #   parTS$alpha <- rep(alpha,nTimes)
  #   parTS$beta <- rep(beta,nTimes)
  } else if (SWGparameterization=='seas'){
    parTS$pdd <- assignSeasPars(parS['pdd.SON'], parS['pdd.DJF'], parS['pdd.MAM'], parS['pdd.JJA'], datInd[["i.ss"]])
    parTS$pwd <- assignSeasPars(parS['pwd.SON'], parS['pwd.DJF'], parS['pwd.MAM'], parS['pwd.JJA'], datInd[["i.ss"]])
    parTS$alpha <- assignSeasPars(parS['alpha.SON'], parS['alpha.DJF'], parS['alpha.MAM'], parS['alpha.JJA'], datInd[["i.ss"]])
    parTS$beta <- assignSeasPars(parS['beta.SON'], parS['beta.DJF'], parS['beta.MAM'], parS['beta.JJA'], datInd[["i.ss"]])
  } else if (SWGparameterization=='har'){
    #Culley 2019 these parameter generators ignore leap years, so for long time series will become out of sync - DM comment: should create new function for applying harmonic to year, dealing with leap years
    parTS$pdd = harmonicFunc(x=seq(1:nTimes),mean=parS['pdd.m'],amp=parS['pdd.amp'],phase.ang = parS['pdd.ang'],k=1,nperiod=365)
    parTS$pwd = harmonicFunc(x=seq(1:nTimes),mean=parS['pwd.m'],amp=parS['pwd.amp'],phase.ang = parS['pwd.ang'],k=1,nperiod=365)
    parTS$alpha = harmonicFunc(x=seq(1:nTimes),mean=parS['alpha.m'],amp=parS['alpha.amp'],phase.ang = parS['alpha.ang'],k=1,nperiod=365)
    parTS$beta = harmonicFunc(x=seq(1:nTimes),mean=parS['beta.m'],amp=parS['beta.amp'],phase.ang = parS['beta.ang'],k=1,nperiod=365)
  }
  
  #Culley 2019 setting 0-1 limits for pdd,pwd.
  parTS$pdd[parTS$pdd>1] = 1.
  parTS$pdd[parTS$pdd<0] = 0.
  
  parTS$pwd[parTS$pwd>1] = 1.
  parTS$pwd[parTS$pwd<0] = 0.
  
  #Culley 2019 non negative limits for alpha,beta
  parTS$alpha[parTS$alpha<.Machine$double.xmin] = .Machine$double.xmin
  parTS$beta[parTS$beta<.Machine$double.xmin] = .Machine$double.xmin
  
  return(parTS)
  
}

#################################

SWGsim.wgen = function(SWGpar,
                          nTimes,
                          randomTerm,
                          auxInfo=NULL){

  # sim occurrence
  simS=Pstatus_WGEN(parPwd=SWGpar$pwd,    # vector of pars for pwd (length = ndays)
                    parPdd=SWGpar$pdd,    # vector of pars for pdd (length = ndays)
                    ndays=nTimes,
                    randomVector = randomTerm$randomVector)

  # sim amounts
  simP=Pamount_WGEN(parAlpha=SWGpar$alpha,         # vector of pars for alpha (length = ndays)
                    parBeta=SWGpar$beta,           # vector of pars for beta (length = ndays)
                    status_ts=simS,            # TS vector of wet/dry statuses-obtained from the output of 'wvar_gen_Pstatus'
                    seed=randomTerm$seed,                       # random seed
                    ndays=nTimes)

  return(simP$sim)
}

Pstatus_WGEN <- function(parPwd,    # vector of pars for pwd (length = nperiod) - The modified cpp code expects vector of length ndays (not nperiod)
                         parPdd,    # vector of pars for pdd (length = nperiod) - The modified cpp code expects vector of length ndays (not nperiod)
                         ndays,
                         randomVector

){

  drywet_TS <- Pstatus_WGEN_cpp(parPwd, parPdd, randomVector, ndays)
  return(drywet_TS)
}

Pamount_WGEN <- function(parAlpha=NULL,         # vector of pars for alpha (length = nperiod) - alpha has to be of length ndays (as per the code)
                         parBeta=NULL,          # vector of pars for beta (length = nperiod)
                         status_ts=NULL,             # TS vector of wet/dry statuses-obtained from the output of 'wvar_gen_Pstatus'
                         seed=NULL,                # random seeds
                         ndays=NULL
){

  set.seed(seed)     # seed seed to fix input needed for rgamma  ---  this creates a challenge for passing in a vector of random numbers
  rain <- vector(mode="numeric", ndays)  #allocate time series
  wet.days<-which(status_ts==1)         #index wet occurance
  rain[wet.days] <- stats::rgamma(parAlpha[wet.days],shape=parAlpha[wet.days],scale=parBeta[wet.days])   #sample rain amount
  syntP <- list(sim=rain,  #
                seed=seed)
  return(syntP)
}

#################################
