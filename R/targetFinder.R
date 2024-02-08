########################################################
## FUNCTION TO SIMULATE SERIES AND DETERMINE FITNESS ###
########################################################

#CONTAINS
  #targetFinder() -
#------------------------------------------------------------------------------------------------
targetFinder<- function(x,               # vector of pars (will change in optim)
                        modelInfo=NULL,
                        modelEnv=NULL,      # tag to link/identify model
                        attSel=NULL,        # attributes selected (vector of strings)
                        attPrim=NULL,       # primary attribute label
                        attInfo=NULL,     # added info regarding attributes (maybe add in attPrim here)!!!!!!!!!!!!!
                        datInd=NULL,
                        randomVector = NULL,
                        randomUnitNormalVector = NULL,
                        target=NULL,        # target locations: desired changes in climate to be simulated, in % relative or abs diff to baseline levels (vector)
                        attObs=NULL,        # observed series attribute values
                        obs=NULL,
                        lambda.mult=NULL,   # lambda multiplier for penalty function
                        simSeed=NULL,
                        wdSeries=NULL,
                        resid_ts=NULL,
                        returnThis = 'objFunc',
                        obj.func = 'SS_absPenalty'
                        #Nw=NULL,            # warmup period in days
                        # N=NULL,             # seeds
                        # seed1=NULL,
                        # seed2=NULL,
                        # seed3=NULL
){

  parS = x

  #SIMULATE SELECTED VARIABLE USING CHOSEN STOCHASTIC MODEL

  sim=switch_simulator(type=modelInfo$simVar,          # what vartype is being simulated
                       parS=parS,
                       modelEnv=modelEnv,
                       randomVector = randomVector,
                       randomUnitNormalVector = randomUnitNormalVector,
                       wdSeries=wdSeries,
                       resid_ts=resid_ts,
                       seed=simSeed,
                       obs=obs)
   
    if(length(which(is.na(sim$sim))) > 0){
    score=-150  #default here
  }else{
    #CALCULATE SELECTED ATTRIBUTE VALUES
    sim.att=attribute.calculator(attSel=attSel,data=sim$sim,datInd=datInd,attInfo=attInfo)

    #RELATING TO BASELINE SERIES
    simPt=unlist(Map(function(type, val,baseVal) simPt.converter.func(type,val,baseVal), attInfo$targetType, as.vector(sim.att),as.vector(attObs)),use.names = FALSE)

    if(returnThis=='sim'){
      return(simPt)
    } else if (returnThis=='resid') {
      weights = rep(1,length(attSel));names(weights) = attSel
      weights[attPrim] = lambda.mult
      resid = weights*(simPt-unlist(target))
      asr = abs(sum(resid))
      if(is.infinite(asr)|is.nan(asr)|is.na(asr)){browser()}
      return(resid)
    } else if (returnThis=='objFunc'){

      #GET OBJECTIVE FUNCTION VALUE ()
      score=objFuncMC(attSel= attSel,     # vector of selected attributes
                      attPrim=attPrim,      # any primary attributes
                      attInfo=attInfo,
                      simPt=simPt,
                      target=target,
                      obj.func=obj.func,   #make this changeable (auto calc lambda)
                      lambda=lambda.mult
      )
    }
  }
  return(score)
}









