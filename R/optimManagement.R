#######################################
##   OPTIMIZER MANAGER FUNCTIONS     ##
#######################################

# CONTAINS
  # gaWrapper() - allows use of suggestions for initial populations
  # screenSuggest() - remove suggestions outside bounds
   # enforceBounds()
      # testBound()
        # outBound()

#-----------------------------------------------

findFixedPars = function(xLo,xHi){
  fixParLoc = fixParVal = fitParLoc = c()
  for (i in 1:length(xLo)){
    if (xLo[i]==xHi[i]){
      fixParLoc = c(fixParLoc,i)
      fixParVal = c(fixParVal,xLo[i])
    } else{
      fitParLoc = c(fitParLoc,i)
    }
  }
return(list(fixParLoc=fixParLoc,fixParVal=fixParVal,fitParLoc=fitParLoc))
}

#-----------------------------------------------

calcParFixedPars = function(x,fixedPars){
  if (!is.null(fixedPars$fixParLoc)){
    xAll = c()
    xAll[fixedPars$fixParLoc] = fixedPars$fixParVal
    xAll[fixedPars$fitParLoc] = x
  } else {
    xAll = x
  }
#  browser()
  return(xAll)
}

#-----------------------------------------------

targetFinderFixPars = function(x,fixedPars=NULL,returnThis='objFunc',...){
    # deal with fixed and fitted pars
  xAll = calcParFixedPars(x,fixedPars)
  target = targetFinder(x=xAll,returnThis=returnThis,...)
  # save obj func value to vector
 if (returnThis=='objFunc'){
   assign("fTrace",c(foreSIGHT_optimizationDiagnosticsEnv$fTrace,
                     target),
          envir = foreSIGHT_optimizationDiagnosticsEnv)
 } else if (returnThis=='resid'){
   assign("fTrace",c(foreSIGHT_optimizationDiagnosticsEnv$fTrace,
                     -sqrt(sum(target^2))),
          envir = foreSIGHT_optimizationDiagnosticsEnv)
 }
  return(target)
}

#-----------------------------------------------

negTargetFinder = function(x,...){
  target = -targetFinder(x=x,...)
  return(target)
}

#-----------------------------------------------

negTargetFinderFixPars = function(x,fixedPars=NULL,...){
  target = -targetFinderFixPars(x=x,fixedPars=fixedPars,...)
  return(target)
}

#-----------------------------------------------
#
foreSIGHT_optimizationDiagnosticsEnv <- new.env(parent = emptyenv())

#-----------------------------------------------

  multiStartOptim = function(optimArgs=NULL,
                             modelInfo=NULL,
                             lambda.mult=NULL,
                             target=NULL,
                             parSuggest=NULL,
                             simSeed=NULL,
                             ...){

  timeStart=Sys.time()

  ##################################

  xLo = modelInfo$minBound
  xHi = modelInfo$maxBound

  fixedPars = findFixedPars(xLo,xHi)

  fBest = 9e9

  if(!is.null(parSuggest)){
    nSugg = nrow(parSuggest)
  } else {
    nSugg = 0
  }

  fMulti = timeMulti = callsMulti = convergedMulti = convergenceCodeMulti = messageMulti = c()
  parsMulti = onBoundsMulti = matrix(nrow=optimArgs$nMultiStart,ncol=length(xLo))
  fTraceMulti = callsTraceMulti = list()

  convergenceCodeSingle = messageSingle = NA

  if(!is.null(optimArgs$seed)){
    seed1=optimArgs$seed
  } else {
    seed1 = 1
  }

  for (r in 1:optimArgs$nMultiStart){

    time1 = Sys.time()
    assign("WG_calls",0,envir = foreSIGHT_optimizationDiagnosticsEnv)
    assign("fTrace",c(),envir = foreSIGHT_optimizationDiagnosticsEnv)
    #foreSIGHT_optimizationDiagnosticsEnv$WG_calls = 0
    #foreSIGHT_optimizationDiagnosticsEnv$fTrace = c()

    print(r)

    seed = seed1 + r - 1

    set.seed(seed) # set the random seed for selecting initial parameter values. note same set of seeds will be used for each target/replicate.

    if (r<=nSugg){
      x0 = parSuggest[r,]
    } else {
      x0 = xLo + stats::runif(length(xLo))*(xHi-xLo)
    }

    if (optimArgs$optimizer=='RGN') {

      outTmp <- RGN::rgn(simFunc=targetFinderFixPars,
                            fixedPars=fixedPars,
                            par = x0[fixedPars$fitParLoc],
                            upper = xHi[fixedPars$fitParLoc],
                            lower = xLo[fixedPars$fitParLoc],
                            simTarget = rep(0.,length(target)),
                            control=optimArgs$RGN.control,
                            lambda.mult=lambda.mult,
                            modelInfo=modelInfo,
                            target=target,
                            returnThis='resid',
                            simSeed=simSeed,
                            ...)

      fSingle = sqrt(2*outTmp$value)
      parsSingle = calcParFixedPars(outTmp$par,fixedPars)
      convergedSingle = (outTmp$convergence)

     } else if (optimArgs$optimizer=='CMAES') {

       set.seed(r)
        outTmp <- cmaes::cma_es(fn=targetFinderFixPars,
                              par = x0[fixedPars$fitParLoc],
                              fixedPars=fixedPars,
                              modelInfo=modelInfo,
                              target=target,
                              lambda.mult=lambda.mult,
                              obj.func=optimArgs$obj.func,
                              simSeed=simSeed,
                         ...,
                         lower = xLo[fixedPars$fitParLoc],
                         upper = xHi[fixedPars$fitParLoc],
                         control=optimArgs$CMAES.control)

        fSingle = -outTmp$value
        if (is.null(outTmp$par)){
          fSingle=9e9
          parSingle = NULL
        } else {
          parsSingle = calcParFixedPars(outTmp$par,fixedPars)
        }

        convergedSingle = (outTmp$convergence==0)
        convergenceCodeSingle = outTmp$convergence

     # } else if (optimArgs$optimizer=='NLSLM') {
     #
     #   outTmp = minpack.lm::nls.lm(par=x0[fixedPars$fitParLoc],
     #                            lower=xLo[fixedPars$fitParLoc],
     #                            upper=xHi[fixedPars$fitParLoc],
     #                            fn=targetFinderFixPars,
     #                            fixedPars=fixedPars,
     #                            control=list(fnscale=1),
     #                            target=target,
     #                            lambda.mult=lambda.mult,
     #                            modelInfo=modelInfo,
     #                            simSeed=simSeed,
     #                            returnThis='resid',
     #                            ...)
     #
     #
     #   fSingle = sqrt(sum(outTmp$fvec^2))
     #
     #   parsSingle = calcParFixedPars(outTmp$par,fixedPars)
     #
     #   convergedSingle = (outTmp$info%in%c(1,2,3,4))
     #   messageSingle = outTmp$message
     #   convergenceCodeSingle = outTmp$info

      } else if (optimArgs$optimizer=='GA') {

      outTmp = GA::ga(type = "real-valued",
                  fitness=targetFinderFixPars,
                  lower = xLo[fixedPars$fitParLoc],
                  upper = xHi[fixedPars$fitParLoc],
                  pcrossover= optimArgs$GA$pcrossover,
                  pmutation=optimArgs$GA.args$pmutation,
                  maxiter=optimArgs$GA.args$maxiter,
                  popSize = optimArgs$GA.args$popSize,
                  maxFitness = optimArgs$GA.args$maxFitness,
                  run=optimArgs$GA.args$run,
                  seed = seed,
                  parallel = optimArgs$GA.args$parallel,
                  keepBest=optimArgs$GA.args$keepBest,
                  suggestions = parSuggest,
                  monitor = FALSE,             #switchback
                  fixedPars=fixedPars,
                  target=target,
                  lambda.mult=lambda.mult,
                  obj.func=optimArgs$obj.func,
                  modelInfo=modelInfo,
                  simSeed=simSeed,
                  ...)

      fSingle = -outTmp@fitnessValue
      parsSingle = calcParFixedPars(outTmp@solution[1,],fixedPars)

      convergedSingle = NA

      } else if (optimArgs$optimizer=='SCE') {

      outTmp = SoilHyP::SCEoptim(
        FUN=targetFinderFixPars,
        par=x0[fixedPars$fitParLoc],
        lower = xLo[fixedPars$fitParLoc],
        upper = xHi[fixedPars$fitParLoc],
        control=optimArgs$SCE.control,
        fixedPars=fixedPars,
        target=target,
        lambda.mult=lambda.mult,
        obj.func=optimArgs$obj.func,
        modelInfo=modelInfo,
        simSeed=simSeed,
        ...)

      fSingle = outTmp$value
      parsSingle = calcParFixedPars(outTmp$par,fixedPars)

      convergedSingle = (outTmp$convergence==0)

    # } else if (optimArgs$optimizer=='BOBYQA') {
    #
    #   outTmp = nloptr::bobyqa(
    #     fn=negTargetFinderFixPars,
    #     x0=x0[fixedPars$fitParLoc],
    #     lower = xLo[fixedPars$fitParLoc],
    #     upper = xHi[fixedPars$fitParLoc],
    #     fixedPars=fixedPars,
    #     target=target,
    #     lambda.mult=lambda.mult,
    #     modelInfo=modelInfo,
    #     simSeed=simSeed,
    #     ...)
    #
    #   fSingle = outTmp$value
    #   parsSingle = calcParFixedPars(outTmp$par,fixedPars)
    #
    #   convergedSingle = (outTmp$convergence>0)
    #   messageSingle = outTmp$message
    #   convergenceCodeSingle = outTmp$convergence

    # } else if (optimArgs$optimizer=='HJK') {
    #
    #   outTmp = dfoptim::hjkb(
    #     fn=targetFinder,
    #     par=x0,
    #     lower = xLo,
    #     upper = xHi,
    #     control=list(maximize=T),
    #     target=target,
    #     lambda.mult=lambda.mult,
    #     modelInfo=modelInfo,
    #     simSeed=simSeed,
    #     ...)
    #
    #   fSingle = -outTmp$value
    #   parsSingle = outTmp$par
    #
    # } else if (optimArgs$optimizer=='NMK') {
    #   outTmp = dfoptim::nmkb(
    #     fn=targetFinderFixPars,
    #     fixedPars=fixedPars,
    #     par=x0[fixedPars$fitParLoc],
    #     lower = xLo[fixedPars$fitParLoc],
    #     upper = xHi[fixedPars$fitParLoc],
    #     control=list(maximize=T),
    #     target=target,
    #     lambda.mult=lambda.mult,
    #     modelInfo=modelInfo,
    #     simSeed=simSeed,
    #     ...)
    #
    #   fSingle = -outTmp$value
    #   parsSingle = calcParFixedPars(outTmp$par,fixedPars)
    #
    #   browser()

    } else if (optimArgs$optimizer=='optim.LBFGSB') {

      outTmp<- stats::optim(par=x0[fixedPars$fitParLoc],
                     fn=targetFinderFixPars,
                     method='L-BFGS-B',
                     lower=xLo[fixedPars$fitParLoc]-1e-6,
                     upper=xHi[fixedPars$fitParLoc]+1e-6,
                     fixedPars=fixedPars,
                     target=target,
                     control=list(fnscale=-1),
                     lambda.mult=lambda.mult,
                     obj.func=optimArgs$obj.func,
                     modelInfo=modelInfo,
                     simSeed=simSeed,
                     ...)

      fSingle = -outTmp$value
      parsSingle = calcParFixedPars(outTmp$par,fixedPars)

      convergedSingle = (outTmp$convergence==0)
      messageSingle = outTmp$message
      convergenceCodeSingle = outTmp$convergence

    # } else if (optimArgs$optimizer=='Rvmmin') {
    #
    #   outTmp<- Rvmmin::Rvmmin(par=x0[fixedPars$fitParLoc],
    #                  fn=targetFinderFixPars,
    #                  lower=xLo[fixedPars$fitParLoc],
    #                  upper=xHi[fixedPars$fitParLoc],
    #                  fixedPars=fixedPars,
    #                  target=target,
    #                  control=list(maximize=T),
    #                  lambda.mult=lambda.mult,
    #                  modelInfo=modelInfo,
    #                  simSeed=simSeed,
    #                  ...)
    #
    #   fSingle = -outTmp$value
    #   parsSingle = calcParFixedPars(outTmp$par,fixedPars)
    #
    #   convergedSingle = (outTmp$convergence==0)
    #   messageSingle = outTmp$message
    #   convergenceCodeSingle = outTmp$convergence
    #
    # } else if (optimArgs$optimizer=='Rcgmin') {
    #
    #   outTmp<- Rcgmin::Rcgmin(par=x0[fixedPars$fitParLoc],
    #                           fn=negTargetFinderFixPars,
    #                           lower=xLo[fixedPars$fitParLoc],
    #                           upper=xHi[fixedPars$fitParLoc],
    #                           fixedPars=fixedPars,
    #                           target=target,
    #                           lambda.mult=lambda.mult,
    #                           modelInfo=modelInfo,
    #                           simSeed=simSeed,
    #                           ...)
    #
    #   fSingle = outTmp$value
    #   parsSingle = calcParFixedPars(outTmp$par,fixedPars)
    #
    #   convergedSingle = (outTmp$convergence==0)
    #   messageSingle = outTmp$message
    #   convergenceCodeSingle = outTmp$convergence

    } else if (optimArgs$optimizer=='NM') {

      outTmp<- dfoptim::nmkb(par=x0[fixedPars$fitParLoc],
                     fn=targetFinderFixPars,
                     lower=xLo[fixedPars$fitParLoc],
                     upper=xHi[fixedPars$fitParLoc],
                     fixedPars=fixedPars,
                     target=target,
                     control=optimArgs$NM.control,
                     lambda.mult=lambda.mult,
                     obj.func=optimArgs$obj.func,
                     modelInfo=modelInfo,
                     simSeed=simSeed,
                     ...)

      fSingle = -outTmp$value
      parsSingle = calcParFixedPars(outTmp$par,fixedPars)

      convergedSingle = (outTmp$convergence==0)
      messageSingle = outTmp$message
      convergenceCodeSingle = outTmp$convergence

    }

    time2 = Sys.time()
    timeSingle=time2-time1    #optimisation runtime

    onBoundsSingle = (abs(parsSingle-xLo)<1e-6) | (abs(parsSingle-xHi)<1e-6)

#    print(parsSingle)

    # calculate best OF value after each function call
    fTrace = foreSIGHT_optimizationDiagnosticsEnv$fTrace
    nTrace = length(fTrace)
    fTraceTmp = c()
    fTraceBest = -999
    callsTraceTmp = 1:nTrace
    for (i in 1:nTrace){
      if ( !is.na(fTrace[i]) & (fTrace[i]>fTraceBest) ){
        fTraceBest = fTrace[i]
      }
      fTraceTmp[i] = fTraceBest
    }
    # remove best OF values when they are repeated
    fTraceTrim = c(fTraceTmp[1])
    callsTraceTrim = c(callsTraceTmp[1])
    for (i in 2:(nTrace-1)){
      if ((fTraceTmp[i]!=fTraceTmp[i-1])|(fTraceTmp[i]!=fTraceTmp[i+1])){
        fTraceTrim = c(fTraceTrim,fTraceTmp[i])
        callsTraceTrim = c(callsTraceTrim,callsTraceTmp[i])
      }
    }

    fMulti[r] = fSingle
    parsMulti[r,] = parsSingle
    timeMulti[r] = timeSingle
    callsMulti[r] = foreSIGHT_optimizationDiagnosticsEnv$WG_calls
    onBoundsMulti[r,] = onBoundsSingle
    convergedMulti[r] = convergedSingle
    convergenceCodeMulti[r] = convergenceCodeSingle
    messageMulti[r] = messageSingle
    if (fSingle<fBest){
      fBest = fSingle
      parsBest = parsSingle
      optOut=outTmp
    }
#    cat(fSingle,fBest,'\n')
    fTraceMulti[[r]] = fTraceTrim
    callsTraceMulti[[r]] = callsTraceTrim
    if(fBest<optimArgs$OFtol){break()}
  }

  timeFin=Sys.time()
  timeRun=timeFin-timeStart    #optimisation runtime

  out=list(par=as.vector(parsBest),
           fitness=as.numeric(fBest),
           seed=simSeed,
           opt=optOut,
           runtime=timeRun,
           fMulti=fMulti,
           parsMulti=parsMulti,
           onBoundsMulti=onBoundsMulti,
           callsMulti=callsMulti,
           timeMulti=timeMulti,
           convergedMulti=convergedMulti,
           convergenceCodeMulti=convergenceCodeMulti,
           messageMulti=messageMulti,
           fTraceMulti=fTraceMulti,
           callsTraceMulti=callsTraceMulti)

  return(out)

}

#-----------------------------------------------
#FUNCTION TO SCREEN SUGGESTED POPULATIONS (screen outside, screen once)
screenSuggest<-function(modelInfo=NULL,
                        modelTag=NULL,
                        parLoc=NULL,        # which pars belong to which model parLoc[[mod]]=c(start,end)
                        suggest=NULL        # suggested pars
){
  nMod=length(modelTag)
  ind=NULL
  if (is.vector(suggest)){suggest=matrix(suggest,nrow=1)}
  for(mod in 1:nMod){
    parSel=suggest[,(parLoc[[mod]][1]:parLoc[[mod]][2])]              #grab par suggestions related to modelTag running
    if (is.vector(parSel)){parSel=matrix(parSel,nrow=1)}
    tmpInd=enforceBounds(suggest=parSel,                              #matrix of suggestions
                         minBound=modelInfo[[modelTag[mod]]]$minBound,
                         maxBound=modelInfo[[modelTag[mod]]]$maxBound)
    #JOIN THE INDICES TOGETHER
    ind=c(ind,tmpInd)
  }

  #REMOVE INDICES DUPLICATES
  ind=unique(ind)

  #REMOVE INAPPROPRIATE SUGGESTIONS
  suggest=suggest[ind,]
  return(suggest)
}

enforceBounds<-function(suggest=NULL, #matrix of suggestions
                        minBound=NULL,
                        maxBound=NULL
){
  nSuggest=nrow(suggest)
  ind=vapply(X=seq(1,nSuggest),FUN=testBound,matPar=suggest,minBound=minBound,maxBound=maxBound,FUN.VALUE=numeric(1))
  ind=ind[which(!is.na(ind))]
  ind
}

testBound<-function(ind=NULL,
                    matPar=NULL,
                    minBound=NULL,
                    maxBound=NULL
){

  out=outBound(ind=ind,inVector=matPar[ind,],minBound=minBound,maxBound=maxBound)
  out
}
# matTest=c(1,2,3,4,4,4,5,5,5,10,11,12);matTest=matrix(matTest,ncol=3,nrow=4,byrow=TRUE)
# minBound=c(1,1,3); maxBound=c(6,8,13)

outBound<-function(ind=NULL,      #index of pars being evaluated
                   inVector=NULL,
                   minBound=NULL,
                   maxBound=NULL
){
  npar=length(minBound)  #No. of pars that should be inside bound
  ntest=length(which((inVector>=minBound)&(inVector<=maxBound)))
  if(ntest==npar){ok=ind}else{ok=NA}
  return(ok)
}

# #
# optpar<- ga(type = "real-valued",
#             fitness=targetFinder,
#             min = modelInfo[[2]]$minBound,
#             max = modelInfo[[2]]$maxBound,
#             pcrossover= gaArgs$pcrossover,
#             pmutation=gaArgs$pmutation,
#             maxiter=gaArgs$maxiter,
#             popSize = gaArgs$popSize,
#             maxFitness = gaArgs$maxFitness,
#             run=gaArgs$run,
#             seed = gaArgs$seed,
#             parallel = gaArgs$parallel,
#             keepBest=gaArgs$keepBest,
#             monitor = FALSE,
#             #BELOW RELATED TO TARGETFINDER()
#             modelTag=modelTag[2],
#             modelInfo=modelInfo[[2]],
#             attSel=attSel[attInd[[2]]],
#             attPrim=attPrim,
#             attInfo=attInfo[[2]],
#             datInd=datInd[[2]],
#             initCalibPars=initCalibPars,
#             target=target[attInd[[2]]],
#             attObs=attObs[attInd[[2]]],
#             lambda.mult=gaArgs$lambda.mult,
#             simSeed=simSeed,
#             wdSeries=wdStatus,
#             resid_ts=resid_ts
# )



#optim=FALSE


# seed1=122
# #.ga.optim.default<-list()
# gaArgs=list(
#   pcrossover= 0.1,
#   pmutation=0.8,
#   maxiter=100,
#   popSize = 100,
#   maxFitness = -0.01,
#   run=300,
#   seed = 1234,
#   parallel = TRUE,
#   keepBest=TRUE
# )

# #TESTER CALL
# optTest=list()
# optTest[[1]]=gaWrapper( gaArgs=gaArgs,  #can specify your own outside
#                       modelTag=modelTag,      # tag to link/identify model
#                       modelInfo=modelInfo,
#                       attSel=attSel,        # attributes selected (vector of strings)
#                       attPrim=attPrim,       # primary attribute label
#                       attInfo=attInfo,       # added info regarding attributes (maybe add in attPrim here)!!!!!!!!!!!!!
#                       datInd=datInd,        # dat ind
#                       initCalibPars=initCalibPars, #vector of pars from initial baseline calibration
#                       target=c(1.05,1.06,1.10),    # target locations: desired changes in climate to be simulated, in % relative or abs diff to baseline levels (vector)
#                       attObs=obs.att,        # observed series attribute values
#                       lambda.mult=2,   # lambda multiplier for penalty function
#                       Nw=100,            # warmup period in days
#                       N=seed1,             # seeds
#                       seed1=seed1,
#                       seed2=seed1,
#                       seed3=seed1
#         )


#Test
# mod=2
# optTest=gaWrapper(gaArgs=.ga.default,          #can specify your own outside
#                     modelTag=modelTag[mod],      # tag to link/identify model
#                     modelInfo=modelInfo[[modelTag[mod]]],
#                     attSel=attSel[attInd[[mod]]],
#                     attPrim=attPrim,
#                     attInfo=attInfo[[modelTag[mod]]],
#                     datInd=datInd[[modelTag[mod]]],
#                     initCalibPars=NULL,
#                     target=targetMat[1,attInd[[2]]],        # target locations: desired changes in climate to be simulated, in % relative or abs diff to baseline levels (vector)
#                     attObs=attObs[attInd[[mod]]],        # observed series attribute values
#                     lambda.mult=1.0,   # lambda multiplier for penalty function
#                     simSeed=1234           # seeds
# )
#
# plot(optTest$opt)
#
# out=switch_simulator(type=modelInfo[[modelTag[mod]]]$simVar,parS=optTest$par,
#                     modelTag=modelTag[mod],modelInfo=modelInfo[[modelTag[mod]]],datInd=datInd[[modelTag[mod]]],
#                     initCalibPars=NULL,wdSeries=NULL,resid_ts=NULL,seed=1234)
#
#   sim.att=attribute.calculator(attSel=attSel[attInd[[mod]]],
#                                    data=out$sim,
#                                    datInd=datInd[[mod]],
#                                    attribute.funcs=attribute.funcs)
