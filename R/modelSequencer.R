###############################
##      MODEL SEQUENCER      ##
###############################

#CONTAINS
  # simulateTarget() - simulates all requested timseries for an individual target location (also outputs attributes and targetSim)

#------------------------------------------------
#FUNCTIONS
simulateTarget<-function(
                    optimArgs=NULL,
                    simVar=NULL,
                    modelTag=NULL,
                    modelInfo=NULL,
                    attSel=NULL,
                    attPrim=NULL,
                    attInfo=NULL,
                    attInd=NULL,
                    datInd=NULL,
                    obs=obs,
                    initCalibPars=NULL,
                    targetLoc=NULL,     #is  a vector  (just 1 target here)
                    attObs=NULL,
                    parLoc=NULL,        #which pars belong to which model parLoc[[mod]]=c(start,end)
                    parSim=NULL,        #pars used to simulate targets so far
                    setSeed=1234,
                    file=NULL,
                    randomUnitNormalVector=NULL,
                    iRepTarg = NULL
                   # resid_ts=NULL    - for other models not currently in play
                    ){

  #LOOP OVER EACH STOCHASTIC MODEL NOMINATED
  nMod=length(modelTag)
  out=list()
  parV=NULL
  objScore=NULL

  #MERGE WITH ANY SUGGESTIONS SUPPLIED IN OPTIMARGS
  if(!is.null(optimArgs$suggestions)){
    parSugg=rbind(parSim,optimArgs$suggestions)
    if(optimArgs$optimizer=='GA'){
      if(dim(parSugg)[1] > optimArgs$popSize){
        parSugg=parSugg[(1:optimArgs$popSize),]  #make sure it doesn't exceed popSize
      }
    }
  } else {
    parSugg=NULL
  }

  #set.seed(setSeed)

  attSim=list()          #Make list to store simulated attributes
  targetSim=list()       #Make list to store simulated attributes(target space converted)
  for(mod in modelTag){
      
    out[[simVar[mod]]] = list()
    
    if (is.null(randomUnitNormalVector)){
      randomVector <- stats::runif(n=datInd[[mod]][[aggNameShort[[obs$timeStep]]]]$nTimes) # Random vector to be passed into weather generator to reduce runtime
    } else {
      randomVector=NULL
    }

     #IF CONDITIONED ON DRY-WET STATUS, populate wdStatus
      switch(simVar[mod], #
             "P" = {wdStatus=NULL},
             "Temp" = {if(modelInfo[[mod]]$WDcondition==TRUE){
                           wdStatus=out[["P"]]$sim
                         }else{
                           wdStatus=NULL
                         }
                       },
                    {wdStatus=NULL}  #default
             )

    # write data to model environment
    #----------------------------------
    write_model_env(envir = foreSIGHT_modelEnv,
                    modelInfo = modelInfo[[mod]],
                    modelTag = mod,
                    datInd = datInd[[mod]]
                    )
    #-----------------------------------

    parMin = modelInfo[[mod]]$minBound
    parMax = modelInfo[[mod]]$maxBound

    if(length(which(parMin==parMax))==length(parMin)){#
      progress(p("    Working on variable ",simVar[mod]),file)
      progress(p("    Parameters specified by user, no optimisation ..."),file)
      
      out[[simVar[mod]]]$sim = simClim(parS=parMin,              
                                    modelTag = mod,
                                    ppTypes=modelInfo$ppTypes,
                                    datInd=datInd[[mod]],
                                    randomTerm = list(randomVector = randomVector,
                                                      randomUnitNormalVector = randomUnitNormalVector,
                                                      seed=setSeed),
                                    obs=obs[[simVar]])
      
      parV=c(parV,parMin)

    }else{
      progress(p("    Working on variable ",simVar[mod]),file)
      progress(p("    Commencing optimisation..."),file)

      #GRAB PAR SUGGESTIONS RELATED TO modelTag
      if(!is.null(parSugg)){
        parSel=parSugg[,(parLoc[[mod]][1]:parLoc[[mod]][2])] #grab par suggestions related to modelTag running
        if (is.vector(parSel)){parSel=matrix(parSel,nrow=1)}
      }else{
        parSel=NULL      #no suggestions to be had
      }

      timeStep = modelInfo[[mod]]$timeStep

      # calculate attribute info prior to start of optimization so it does not need to be recalculated each call to target finder
      
      
      attInfo$attCalcInfo = list()
      
      aggs = unique(attInfo$aggType)

      for (i in aggs){
        attInfo[[mod]]$attCalcInfo = attribute.calculator.setup(attCalcInfo=attInfo[[mod]]$attCalcInfo,
                                                         attSel=attSel[attInd[[simVar[mod]]]],
                                                         datInd=datInd[[mod]][[i]])
      }

      optTest = multiStartOptim(optimArgs=optimArgs,
                                modelEnv = foreSIGHT_modelEnv,
                                modelInfo=modelInfo[[mod]],
                                modelTag=mod,
                                attSel=attSel[attInd[[simVar[mod]]]],
                                attPrim=attPrim,
                                attInfo=attInfo[[mod]],
                                datInd=datInd[[mod]],
                                randomVector = randomVector,
                                randomUnitNormalVector = randomUnitNormalVector,
                                parSuggest=parSel,
                                target=targetLoc[attInd[[simVar[mod]]]],
                                attObs=attObs[attInd[[simVar[mod]]]],
                                obs=obs,
                                lambda.mult=optimArgs$lambda.mult,
                                simSeed=setSeed,
                                iRepTarg = iRepTarg,
                                wdSeries=wdStatus,   #selecting rainfall  if needed
                                resid_ts=NULL)

      out[[simVar[mod]]] = append(out[[simVar[mod]]],optTest)

      #progress(p("    Best fitness: ",signif(optTest$fitness,digits=5), ". Optimisation stopped at iter ",nIter),file)
      #progress(p("    Note:",signif(summary(optTest$opt)$fitness,digits=5)),file)

      out[[simVar[mod]]]$sim=simClim(parS=optTest$par,              #RAIN SELECTED
                                     modelTag = mod,
                                     ppTypes=modelInfo$ppTypes,
                                     datInd=datInd[[mod]],
                                     randomTerm = list(randomVector = randomVector,
                                                       randomUnitNormalVector = randomUnitNormalVector,
                                                       seed=optTest$seed),
                                     obs=obs)
      parV=c(parV,optTest$par)

    }
    
    aggPeriods = names(datInd[[mod]])
    nagg = length(aggPeriods)
    timeStep = obs$timeStep 

    
      data = list(times=datInd[[mod]][[aggNameShort[[timeStep]]]]$times,
                      timeStep=timeStep)
      data[[simVar[mod]]] = out[[simVar[mod]]]$sim

      sim.att = aggregate_calculate_attributes(varList=simVar[mod],
                                               aggList=aggPeriods,
                                               data=data,
                                               attSel=attSel[attInd[[simVar[mod]]]],
                                               datInd=datInd[[mod]],
                                               attInfo=attInfo[[mod]])
      attSim[[mod]]=sim.att        #store simulated attributes in list
      
      simPt=unlist(Map(function(type, val,baseVal) simPt.converter.func(type,val,baseVal), attInfo[[mod]]$targetType, sim.att,attObs[attInd[[simVar[mod]]]]),use.names = FALSE)
      names(simPt)=attSel[attInd[[simVar[mod]]]]
      targetSim[[mod]]=simPt             #Store in list

      # dist=eucDist(target=targetLoc[attInd[[mod]]],simPt=simPt)
      # progress(paste0("    Euc Dist ",signif(dist,4)),file)
      #
      # primInd=which(attInfo[[modelTag[mod]]]$primType==TRUE)
      # penalty.score=penaltyFunc_basic(target=targetLoc[attInd[[mod]]][primInd],simPt=simPt[primInd],lambda=optimArgs$lambda.mult[attInfo[[modelTag[mod]]]$primMult])
      # progress(paste0("    Penalty ",signif(penalty.score,4)),file)
      #
      # progress(paste("    target - ",paste(attPrim,": ",signif(targetLoc[attInd[[mod]]][primInd],digits=4),collapse = ", ",sep=""),sep=''),file)
      # progress(paste("    simpt - ",paste(attPrim,": ",signif(simPt[attInd[[mod]]][primInd],digits=4),collapse = ", ",sep=""),sep=''),file)
      # progress(paste("    lambda - ",paste(attPrim,": ",signif(optimArgs$lambda.mult[attInfo[[modelTag[mod]]]$primMult],digits=4),collapse = ", ",sep=""),sep=''),file)

      score=objFuncMC(attSel= attSel[attInd[[simVar[mod]]]],     # vector of selected attributes
                      attPrim=attPrim,      # any primary attributes
                      attInfo=attInfo,
                      simPt=simPt,
                      target=targetLoc[attInd[[simVar[mod]]]],
                      obj.func=optimArgs$obj.func,
                      lambda=optimArgs$lambda.mult)

      #CONFIRMING SCORE FOR SIM SERIES
      progress(paste0("    Variable ",simVar[mod]," final sim series fitness: ",signif(score,4)),file)

      objScore=c(objScore,score)

      #CALCULATE SIM ATTRIBUTES HERE (ABSOLUTE AND TARGET SPACE)
      out[[simVar[mod]]]$parS=parV
      out[[simVar[mod]]]$score=objScore
  }  #end model loop

  names(attSim) = names(targetSim) = NULL
  out$attSim=unlist(attSim)[attSel]        # unlist,relist & make sure order is correct
  progress(paste("    Attributes Simulated - ",paste(attSel,": ",signif(out$attSim,digits=4),collapse = ", ",sep=""),sep=''),file)

  out$targetSim=unlist(targetSim)[attSel]  # unlist,relist & make sure order is correct
  progress(paste("    Target Simulated - ",paste(attSel,": ",signif(out$targetSim,digits=4),collapse = ", ",sep=''),sep=""),file)
  
  

  return(out)
}


#TESTER
# tmp=simulateTarget(optimArgs=optimArgs,
#                    simVar=simVar,
#                    modelTag=modelTag,
#                    modelInfo=modelInfo,
#                    attSel=attSel,
#                    attPrim=attPrim,
#                    attInfo=attInfo,
#                    attInd=attInd,
#                    datInd=datInd,
#                    initCalibPars=NULL,
#                    targetLoc=targetMat[1,],     #is  a vector  (just 1 target here)
#                    attObs=attObs,
#                    lambda.mult=1.0,
#                    setSeed=1234)
#
# tmp[[simVar[mod]]]
# tmp$attSim
