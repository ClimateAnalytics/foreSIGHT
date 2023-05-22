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

targetFinderFixPars = function(x,fixedPars=NULL,...){
  # deal with fixed and fitted pars
  xAll = calcParFixedPars(x,fixedPars)
  target = targetFinder(x=xAll,...)
#  print(x)
#  print(target)
#  browser()
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

# multiStartOptim = function(optimArgs=NULL,
#                            modelEnv = NULL,
#                            modelInfo=NULL,     # information related to modelTags
#                            attSel=NULL,        # attributes selected (vector of strings)
#                            attPrim=NULL,       # primary attribute label
#                            attInfo=NULL,       # added info regarding attributes
#                            datInd=NULL,
#                            randomVector = NULL,
#                            randomUnitNormalVector = NULL,
#                            parSuggest=NULL,    # paramater suggestions
#                            target=NULL,        # target locations: desired changes in climate to be simulated, in % relative or abs diff to baseline levels (vector)
#                            attObs=NULL,        # observed series attribute values
#                            lambda.mult=NULL,   # lambda multiplier for penalty function
#                            simSeed=NULL,       # seeds
#                            wdSeries=NULL,
#                            resid_ts=NULL
#                            ){

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

  ######## THESE DEFAULT SETTINGS SHOULD BE INCORPORATED INTO optimArgsdefault() in default_parameters.R . This will require a bit of work to separate GA from RGN settings


  if(!is.null(optimArgs$iterMax)){
    iterMax = optimArgs$iterMax
  } else {
    iterMax = 100
  }

  if(!is.null(optimArgs$suggestions)){
    sugg = optimArgs$suggestions
    nSugg = nrow(sugg)
  } else {
    nSugg = 0
  }

  fMulti = timeMulti = callsMulti = c()
  parsMulti = onBoundsMulti = matrix(nrow=optimArgs$nMultiStart,ncol=length(xLo))

  for (r in 1:optimArgs$nMultiStart){

    time1 = Sys.time()
    assign("WG_calls",0,envir = .GlobalEnv)

    print(r)

    if (optimArgs$optimizer!='GA') {
      if (r<=nSugg){
        x0 = sugg[r,]
      } else {
        set.seed(r) # set the random seed for selecting initial parameter values. note same set of seeds will be used for each target/replicate.
        x0 = xLo + runif(length(xLo))*(xHi-xLo)
      }
    } else {
      if(!is.null(optimArgs$seed)){
        seed=optimArgs$seed[r]
      } else {
        seed = r
      }
    }

    if (optimArgs$optimizer=='RGN') {

#      browser()

      # outTmp <- rgn(simFunc=targetFinderFixPars,
      #                       fixedPars=fixedPars,
      #                       x0 = x0[fixedPars$fitParLoc],
      #                       xHi = xHi[fixedPars$fitParLoc],
      #                       xLo = xLo[fixedPars$fitParLoc],
      #                       simTarget = unlist(target),
      #                       weights=1.+lambda.mult,
      #                       modelInfo=modelInfo,
      #                       target=target,
      #                       returnThis='sim',
      #                       simSeed=simSeed,
      #                       ...)

      outTmp <- rgn(simFunc=targetFinderFixPars,
                            fixedPars=fixedPars,
                            x0 = x0[fixedPars$fitParLoc],
                            xHi = xHi[fixedPars$fitParLoc],
                            xLo = xLo[fixedPars$fitParLoc],
                            simTarget = rep(0.,length(target)),
                            lambda.mult=lambda.mult,
                            modelInfo=modelInfo,
                            target=target,
                            returnThis='resid',
                            simSeed=simSeed,
                            ...)

      fSingle = sqrt(2*outTmp$info$f)
      parsSingle = calcParFixedPars(outTmp$x,fixedPars)

     } else if (optimArgs$optimizer=='CMAES') {

#       browser()

       set.seed(r)
        outTmp <- cmaes::cma_es(fn=targetFinderFixPars,
                              par = x0[fixedPars$fitParLoc],
                              fixedPars=fixedPars,
                              modelInfo=modelInfo,
                              target=target,
                              lambda.mult=lambda.mult,
                              simSeed=simSeed,
                         ...,
                         lower = xLo[fixedPars$fitParLoc],
                         upper = xHi[fixedPars$fitParLoc],
                         control=list(fnscale=-1,keep.best=T,stopfitness=1e-5))#,maxit=100))

       # outTmp1 <- cmaesr::cmaes(objective.fun=targetFinderFixPars,
       #                         start.point = x0[fixedPars$fitParLoc],
       #                         fixedPars=fixedPars,
       #                         modelInfo=modelInfo,
       #                         target=target,
       #                         lambda.mult=lambda.mult,
       #                         simSeed=simSeed,
       #                         ...)


#        browser()

        fSingle = -outTmp$value
        if (is.null(outTmp$par)){
          fSingle=9e9
          parSingle = NULL
        } else {
          parsSingle = calcParFixedPars(outTmp$par,fixedPars)
        }

#        browser()

     } else if (optimArgs$optimizer=='NLSLM') {

       outTmp = minpack.lm::nls.lm(par=x0[fixedPars$fitParLoc],
                                lower=xLo[fixedPars$fitParLoc],
                                upper=xHi[fixedPars$fitParLoc],
                                fn=targetFinderFixPars,
                                fixedPars=fixedPars,
                                control=list(fnscale=1),
                                target=target,
                                lambda.mult=lambda.mult,
                                modelInfo=modelInfo,
                                simSeed=simSeed,
                                returnThis='resid',
                                ...)


       fSingle = sqrt(sum(outTmp$fvec^2))

       parsSingle = calcParFixedPars(outTmp$par,fixedPars)

      } else if (optimArgs$optimizer=='GA') {

      # outTmp = ga(type = "real-valued",
      #             fitness=targetFinder,
      #             lower = xLo,
      #             upper = xHi,
      #             pcrossover= optimArgs$pcrossover,
      #             pmutation=optimArgs$pmutation,
      #             maxiter=optimArgs$maxiter,
      #             popSize = optimArgs$popSize,
      #             maxFitness = optimArgs$maxFitness,
      #             run=optimArgs$run,
      #             seed = r,
      #             parallel = optimArgs$parallel,
      #             keepBest=optimArgs$keepBest,
      #             suggestions = parSuggest,
      #             monitor = FALSE,             #switchback
      #             target=target,
      #             lambda.mult=lambda.mult,
      #             modelInfo=modelInfo,
      #             simSeed=simSeed,
      #             ...)
      #
      # fSingle = -outTmp@fitnessValue
      # parsSingle = outTmp@solution[1,]

      outTmp = ga(type = "real-valued",
                  fitness=targetFinderFixPars,
                  lower = xLo[fixedPars$fitParLoc],
                  upper = xHi[fixedPars$fitParLoc],
                  pcrossover= optimArgs$pcrossover,
                  pmutation=optimArgs$pmutation,
                  maxiter=optimArgs$maxiter,
                  popSize = optimArgs$popSize,
                  maxFitness = optimArgs$maxFitness,
                  run=optimArgs$run,
                  seed = r,
                  parallel = optimArgs$parallel,
                  keepBest=optimArgs$keepBest,
                  suggestions = parSuggest,
                  monitor = FALSE,             #switchback
                  fixedPars=fixedPars,
                  target=target,
                  lambda.mult=lambda.mult,
                  modelInfo=modelInfo,
                  simSeed=simSeed,
                  ...)

      fSingle = -outTmp@fitnessValue
      parsSingle = calcParFixedPars(outTmp@solution[1,],fixedPars)

      } else if (optimArgs$optimizer=='SCE') {

      # outTmp = SoilHyP::SCEoptim(
      #             FUN=targetFinder,
      #             par=x0,
      #             lower = xLo,
      #             upper = xHi,
      #             control=list(fnscale=-1,initsample='random'),
      #             target=target,
      #             lambda.mult=lambda.mult,
      #             modelInfo=modelInfo,
      #             simSeed=simSeed,
      #             ...)
      #
      # fSingle = outTmp$value
      # parsSingle = outTmp$par

      outTmp = SoilHyP::SCEoptim(
        FUN=targetFinderFixPars,
        par=x0[fixedPars$fitParLoc],
        lower = xLo[fixedPars$fitParLoc],
        upper = xHi[fixedPars$fitParLoc],
        control=list(fnscale=-1,initsample='random'),
        fixedPars=fixedPars,
        target=target,
        lambda.mult=lambda.mult,
        modelInfo=modelInfo,
        simSeed=simSeed,
        ...)

      fSingle = outTmp$value
      parsSingle = calcParFixedPars(outTmp$par,fixedPars)

    } else if (optimArgs$optimizer=='BOBYQA') {

      # outTmp = nloptr::bobyqa(
      #   fn=negTargetFinder,
      #   x0=x0,
      #   lower = xLo,
      #   upper = xHi,
      #   target=target,
      #   lambda.mult=lambda.mult,
      #   modelInfo=modelInfo,
      #   simSeed=simSeed,
      #   ...)
      #
      # fSingle = outTmp$value
      # parsSingle = outTmp$par

      outTmp = nloptr::bobyqa(
        fn=negTargetFinderFixPars,
        x0=x0[fixedPars$fitParLoc],
        lower = xLo[fixedPars$fitParLoc],
        upper = xHi[fixedPars$fitParLoc],
        fixedPars=fixedPars,
        target=target,
        lambda.mult=lambda.mult,
        modelInfo=modelInfo,
        simSeed=simSeed,
        ...)

      fSingle = outTmp$value
      parsSingle = calcParFixedPars(outTmp$par,fixedPars)

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

      # outTmp<- optim(par=x0,
      #                fn=targetFinder,
      #                method='L-BFGS-B',
      #                lower=xLo-1e-6,
      #                upper=xHi+1e-6,
      #                target=target,
      #                control=list(fnscale=-1),
      #                lambda.mult=lambda.mult,
      #                modelInfo=modelInfo,
      #                simSeed=simSeed,
      #                ...)
      #
      # fSingle = -outTmp$value
      # parsSingle = outTmp$par

      outTmp<- optim(par=x0[fixedPars$fitParLoc],
                     fn=targetFinderFixPars,
                     method='L-BFGS-B',
                     lower=xLo[fixedPars$fitParLoc]-1e-6,
                     upper=xHi[fixedPars$fitParLoc]+1e-6,
                     fixedPars=fixedPars,
                     target=target,
                     control=list(fnscale=-1),
                     lambda.mult=lambda.mult,
                     modelInfo=modelInfo,
                     simSeed=simSeed,
                     ...)

      fSingle = -outTmp$value
      parsSingle = calcParFixedPars(outTmp$par,fixedPars)

#      browser()

    # } else if (optimArgs$optimizer=='optim.NM') {
    #
    #   outTmp<- optim(par=x0,
    #                  fn=targetFinder,
    #                  method='Nelder-Mead',
    #                  target=target,
    #                  control=list(fnscale=-1),
    #                  lambda.mult=lambda.mult,
    #                  modelInfo=modelInfo,
    #                  simSeed=simSeed,
    #                  ...)
    #
    #   fSingle = -outTmp$value
    #   parsSingle = outTmp$par

    }

    time2 = Sys.time()
    timeSingle=time2-time1    #optimisation runtime

    onBoundsSingle = (abs(parsSingle-xLo)<1e-6) | (abs(parsSingle-xHi)<1e-6)

    print(parsSingle)

#    browser()

    fMulti[r] = fSingle
    parsMulti[r,] = parsSingle
    timeMulti[r] = timeSingle
    callsMulti[r] = WG_calls
    onBoundsMulti[r,] = onBoundsSingle
    if (fSingle<fBest){
      fBest = fSingle
      parsBest = parsSingle
      optOut=outTmp
    }

    cat(fSingle,fBest,'\n')

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
           timeMulti=timeMulti)

  return(out)

}

#-----------------------------------------------

gaWrapper<-function(gaArgs=NULL,        # can specify your own outside
                    modelEnv=NULL,
                    modelInfo=NULL,     # information related to modelTags
                    attSel=NULL,        # attributes selected (vector of strings)
                    attPrim=NULL,       # primary attribute label
                    attInfo=NULL,       # added info regarding attributes
                    datInd=NULL,
                    randomVector = NULL,
                    randomUnitNormalVector = NULL,
                    parSuggest=NULL,    # paramater suggestions
                    target=NULL,        # target locations: desired changes in climate to be simulated, in % relative or abs diff to baseline levels (vector)
                    attObs=NULL,        # observed series attribute values
                    lambda.mult=NULL,   # lambda multiplier for penalty function
                    simSeed=NULL,       # seeds
                    wdSeries=NULL,
                    resid_ts=NULL
                    ){

  #MEMOISE YAY OR NAY - TBC
  #USER SPECIFIED PARALLEL CONTROLS
  timeStart=Sys.time()

  if(!is.null(gaArgs$nMultiStart)){
    nMultiStart = gaArgs$nMultiStart
  } else {
    nMultiStart = 1
  }

  if(!is.null(gaArgs$seed)){
    seedList=gaArgs$seed
  } else {
    seedList = 1:nMultiStart
  }

  fBest = -9e9

  for (r in 1:nMultiStart){     ####### NOTE multi-starts should be done in separate wrapper

    print(r)
    print(seedList[r])

    optparTmp<- ga(type = "real-valued",
                fitness=targetFinder,
                lower = modelInfo$minBound,
                upper = modelInfo$maxBound,
                pcrossover= gaArgs$pcrossover,
                pmutation=gaArgs$pmutation,
                maxiter=gaArgs$maxiter,
                popSize = gaArgs$popSize,
                maxFitness = gaArgs$maxFitness,
                run=gaArgs$run,
                seed = seedList[r],
                parallel = gaArgs$parallel,
                keepBest=gaArgs$keepBest,
                suggestions = parSuggest,
                monitor = FALSE,             #switchback
                #BELOW RELATED TO TARGETFINDER()
                modelInfo=modelInfo,
                modelEnv=modelEnv,
                attSel=attSel,
                attPrim=attPrim,
                attInfo=attInfo,
                datInd=datInd,
                randomVector = randomVector,
                randomUnitNormalVector = randomUnitNormalVector,
                target=target,
                attObs=attObs,
                lambda.mult=lambda.mult,
                simSeed=simSeed,
                wdSeries=wdSeries,
                resid_ts=resid_ts
    )

    if (optparTmp@fitnessValue>fBest){
      fBest = optparTmp@fitnessValue
      optpar=optparTmp
    }

    print(paste(optparTmp@fitnessValue,fBest))

  }

  timeFin=Sys.time()
  timeRun=timeFin-timeStart    #optimisation runtime
  #print(summary(optpar)$fitness)

  out=list(par=as.vector(optpar@solution[1,]),
           fitness=as.numeric(optpar@fitnessValue),
           seed=simSeed,
           opt=optpar,
           runtime=timeRun)

  return(out)

}

#-----------------------------------------------

rgnWrapper<-function(rgnArgs=NULL,        # can specify your own outside
                    modelEnv=NULL,
                    modelInfo=NULL,     # information related to modelTags
                    attSel=NULL,        # attributes selected (vector of strings)
                    attPrim=NULL,       # primary attribute label
                    attInfo=NULL,       # added info regarding attributes
                    datInd=NULL,
                    randomVector = NULL,
                    randomUnitNormalVector = NULL,
                    parSuggest=NULL,    # parameter suggestions
                    target=NULL,        # target locations: desired changes in climate to be simulated, in % relative or abs diff to baseline levels (vector)
                    attObs=NULL,        # observed series attribute values
                    lambda.mult=NULL,   # lambda multiplier for penalty function
                    simSeed=NULL,       # seeds
                    wdSeries=NULL,
                    resid_ts=NULL){

  timeStart=Sys.time()

  ##################################

  xLo = modelInfo$minBound
  xHi = modelInfo$maxBound

  fBest = 9e9

  ######## THESE DEFAULT SETTINGS SHOULD BE INCORPORATED INTO optimArgsdefault() in default_parameters.R . This will require a bit of work to separate GA from RGN settings


  if(!is.null(rgnArgs$iterMax)){
    iterMax = rgnArgs$iterMax
  } else {
    iterMax = 100
  }

  if(!is.null(rgnArgs$rgnSettings$suggestions)){
    sugg = rgnArgs$rgnSettings$suggestions
    nSugg = nrow(sugg)
  } else {
    nSugg = 0
  }

  fMulti = c()
  for (r in 1:rgnArgs$nMultiStart){

    print(r)
    if (r<=nSugg){
      x0 = sugg[r,]
    } else {
      set.seed(r) # set the random seed for selecting initial parameter values. note same set of seeds will be used for each target/replicate.
      x0 = xLo + runif(length(xLo))*(xHi-xLo)
    }
  rgnOutTmp <- rgn_fixPars(simFunc=targetFinder,
                     x0 = x0,
                     xHi = xHi,
                     xLo = xLo,
                     simTarget = unlist(target),
                     weights=lambda.mult,
                     info =rgnInfoType,
                     cnv = setDefaultRgnConvergeSettings(dump=0, fail=0,iterMax = iterMax),
                     modelInfo=modelInfo,
                     modelEnv=modelEnv,
                     attSel=attSel,
                     attPrim=attPrim,
                     attInfo=attInfo,
                     datInd=datInd,
                     randomVector = randomVector,
                     randomUnitNormalVector = randomUnitNormalVector,
                     simSeed=simSeed,
                     wdSeries=wdSeries,
                     resid_ts=resid_ts,
                     attObs=attObs,target=target,return_sim_only=T)

    fMulti[r] = rgnOutTmp$info$f
    if (rgnOutTmp$info$f<fBest){
      fBest = rgnOutTmp$info$f
      rgnOut=rgnOutTmp
    }

    cat(rgnOutTmp$info$f,fBest,'\n')

  }

  timeFin=Sys.time()
  timeRun=timeFin-timeStart    #optimisation runtime

  out=list(par=as.vector(rgnOut$x),
           fitness=as.numeric(rgnOut$info$f),
           seed=simSeed,
           opt=rgnOut,
           runtime=timeRun,
           fMulti=fMulti)

  return(out)

}

#FUNCTION TO SCREEN SUGGESTED POPULATIONS (screen outside, screen once)
screenSuggest<-function(modelInfo=NULL,
                        modelTag=NULL,
                        parLoc=NULL,        # which pars belong to which model parLoc[[mod]]=c(start,end)
                        suggest=NULL        # suggested pars
){
  nMod=length(modelTag)
  ind=NULL
  for(mod in 1:nMod){
    parSel=suggest[,(parLoc[[mod]][1]:parLoc[[mod]][2])]              #grab par suggestions related to modelTag running
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
  ntest=length(which((inVector>minBound)&(inVector<maxBound)))
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
