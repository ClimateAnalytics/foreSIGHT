######################################################

combine_sims = function(sim.1,var.1,sim.2,var.2,combination.reps=T){
  
  sim.new = list()
  
  if(identical(sim.1$simDates,sim.2$simDates)){
    sim.new$simDates = sim.1$simDates
  } else {
    print('simDates differ')
    return()
  }
  
  sim.new$controlFile.1 = sim.1$controlFile
  sim.new$controlFile.2 = sim.2$controlFile
  
  expSpace.1 = sim.1$expSpace
  expSpace.2 = sim.2$expSpace
  expSpace.new = list()
  # if (expSpace.1$targetType==expSpace.PET$targetType){
  #   expSpace.P$targetType==expSpace.PET$targetType
  # }
  expSpace.new$attPerturb = c(expSpace.1$attPerturb,expSpace.2$attPerturb)
  expSpace.new$attPerturbSamp = c(expSpace.1$attPerturbSamp,expSpace.2$attPerturbSamp)
  expSpace.new$attPerturbMin = c(expSpace.1$attPerturbMin,expSpace.2$attPerturbMin)
  expSpace.new$attPerturbMax = c(expSpace.1$attPerturbMax,expSpace.2$attPerturbMax)
  expSpace.new$targetType = c(expSpace.1$targetType,expSpace.2$targetType)
  
  tar.1 = 1:nrow(expSpace.1$targetMat)
  tar.2 = 1:nrow(expSpace.2$targetMat)
  
  atts.1 = colnames(expSpace.1$targetMat)
  atts.2 = colnames(expSpace.2$targetMat)
  
  atts.total = c(atts.1,atts.2)
  
  common = intersect(atts.1,atts.2)
  
  if(length(common)>0){
    print('sim1 and sim2 cannot have common attributes in expSpaces')
    return()
  }
  
  nTar.total = length(tar.1)*length(tar.2)
  targetMat = matrix(nrow=nTar.total,ncol=length(atts.total))
  tar.index.1 = tar.index.2 = rep(NA,nrow(targetMat))
  tar = 0
  for (t.2 in tar.2){
    for (t.1 in tar.1){
      tar = tar + 1
      tar.index.1[tar] = t.1
      tar.index.2[tar] = t.2
      targetMat[tar,] = c(unlist(expSpace.1$targetMat[t.1,]),
                          unlist(expSpace.2$targetMat[t.2,]))
    }
  }
  targetMat = data.frame(targetMat)
  colnames(targetMat) = atts.total
  expSpace.new$targetMat = targetMat
  sim.new$expSpace = expSpace.new
  
  nReps.1 = length(which(grepl('Rep',names(sim.1))))
  nReps.2 = length(which(grepl('Rep',names(sim.2))))
  if (combination.reps){
    nReps.total = nReps.1*nReps.2
    rep.index.1 = rep.index.2 = rep(NA,nReps.total)
    rep = 0
    for (r.2 in 1:nReps.2){
      for (r.1 in 1:nReps.1){
        rep = rep + 1
        rep.index.1[rep] = r.1
        rep.index.2[rep] = r.2
      }
    }
  } else {
    if (!nReps.1==nReps.2){
      print('require nReps1==nReps2')
      return()
    }
    nReps.total = nReps.1
    rep.index.1 = rep.index.2 = rep(NA,nReps.total)
    for (r in 1:nReps.1){
      rep.index.1[r] = r
      rep.index.2[r] = r
    }
  }
  
  
  for (r in 1:nReps.total){
    repName = paste0('Rep',r)
    repName.1 = paste0('Rep',rep.index.1[r])
    repName.2 = paste0('Rep',rep.index.2[r])
    sim.new[[repName]] = list()
    for (t in 1:nTar.total){
      tarName = paste0('Target',t)
      tarName.1 = paste0('Target',tar.index.1[t])
      tarName.2 = paste0('Target',tar.index.2[t])
      sim.new[[repName]][[tarName]][[var.1]] = sim.1[[repName.1]][[tarName.1]][[var.1]]
      sim.new[[repName]][[tarName]][[var.2]] = sim.2[[repName.2]][[tarName.2]][[var.2]]
    }
  }
  
  return(sim.new)
  
}

######################################################

add_obs_var_to_sim = function(sim,var,data){

  sim.new = sim
  
  nReps = length(which(grepl('Rep',names(sim))))
  nTar = length(sim[[1]])
  
  for (r in 1:nReps){
    repName = paste0('Rep',r)
    for (t in 1:nTar){
      tarName = paste0('Target',t)
      sim.new[[repName]][[tarName]][[var]] = list(sim=data)      
    }
  }

  sim.new$mergeInfo = paste0('added ',var)
  
  return(sim.new)  

}


######################################################

shuffle_sim = function(sim,clim,attPerturb='P_day_all_tot_dwellTime',
                       targetVals,
                       targetType='frac',seed=1,
                       annAR1coeffList = seq(-0.2,0.9,0.01),
                       cSel='mean'){
  
  # browser()
  
  perturb.varname = get.attribute.varType(attPerturb)
  
  varNames = unlist(lapply(colnames(sim$expSpace$targetMat),get.attribute.varType))
  varNames = unique(unlist(strsplit(varNames,'/')))
  other.varnames = varNames[varNames!=perturb.varname]
  
  expSpace = sim$expSpace
  expSpace.new = list()
  expSpace.new$attPerturb = c(expSpace$attPerturb,attPerturb)
  expSpace.new$attHold = expSpace$attHold
  expSpace.new$targetType = c(expSpace$targetType,targetType)
  expSpace.new$attTied = expSpace$attTied
  
  targetMat.orig = expSpace$targetMat
  
  tar.1 = 1:nrow(targetMat.orig)
  tar.2 = 1:length(targetVals)
  nTar.total = length(tar.1)*length(tar.2)
  
  atts.1 = colnames(targetMat.orig)
  atts.total = c(atts.1,attPerturb)
  nAtt.Total = length(atts.total)
  
  targetMat = matrix(nrow=nTar.total,ncol=nAtt.Total)
  tar.index.1 = tar.index.2 = rep(NA,nrow(targetMat))
  tar = 0
  for (t.2 in tar.2){
    for (t.1 in tar.1){
      tar = tar + 1
      tar.index.1[tar] = t.1
      tar.index.2[tar] = t.2
      targetMat[tar,] = c(unlist(targetMat.orig[t.1,]),targetVals[t.2])
    }
  }
  targetMat = data.frame(targetMat)
  colnames(targetMat) = atts.total
  expSpace.new$targetMat = targetMat
  
  nReps = length(which(grepl('Rep',names(sim))))
  
  d=dim(clim_ref[[perturb.varname]])
  multisite = FALSE
  if (!is.null(d)){
    if (d>1){
      multisite = TRUE
    }
  }
  
  if (multisite){
    if (!is.null(cSel)){
      if (cSel=='mean'){
        clim_ref[[perturb.varname]] = apply(clim_ref[[perturb.varname]],1,mean)
      } else {
        clim_ref[[perturb.varname]] = clim_ref[[perturb.varname]][,cSel]
      }      
    }
  }
    
  att.clim = calculateAttributes(clim_ref,attPerturb)
  
  # browser()
  
  att.sim.multi = list()
  sim.new = list()
  for (r in 1:nReps){
    repName = paste0('Rep',r)
    print(repName)
    sim.new[[repName]] = list()
    att.sim.multi[[repName]] = list()
    for (t in 1:nTar.total){
      tarName = paste0('Target',t)
      print(tarName)
      tarName.1 = paste0('Target',tar.index.1[t])
      if(is.null(att.sim.multi[[repName]][[tarName.1]])){
        att.sim.multi[[repName]][[tarName.1]] = list()
        var.sim = sim[[repName]][[tarName.1]][[perturb.varname]]$sim
        times = sim$simDates
        clim_sim = list(times=times)
        clim_sim[[perturb.varname]]=var.sim
        if (multisite){
          if (!is.null(cSel)){
            if (cSel=='mean'){
              clim_sim[[perturb.varname]] = apply(clim_sim[[perturb.varname]],1,mean)
            } else {
              clim_sim[[perturb.varname]] = clim_sim[[perturb.varname]][,cSel]
            }      
          }
        }
        att.sim.multi[[repName]][[tarName.1]] = c()
        for (i in 1:length(annAR1coeffList)){
          annAR1coeff = annAR1coeffList[i]
          var.sim.new = pp.annShuffle(var.sim,times,annAR1coeff,seed)
          clim_sim = list(times=times); clim_sim[[perturb.varname]]=var.sim.new
          att.sim.multi[[repName]][[tarName.1]][i] = calculateAttributes(clim_sim,attPerturb)
        }      
      } 
      
      # sim.new[[repName]][[tarName]][[perturb.varname]] = sim[[repName]][[tarName.1]][[perturb.varname]]

      att.target = targetMat[t,attPerturb]*att.clim
      print(att.target)
      
      abs_diff = abs(att.sim.multi[[repName]][[tarName.1]] - att.target)
      i = min(which(abs_diff == min(abs_diff)))

      #med <- median(i)  # Interpolated median
      # Find the actual value in x closest to it
      #i <- i[which.min(abs(i - med))]
      
      # par(mfrow=c(1,1))
      # ylim = c(min(att.sim.multi[[repName]][[tarName.1]],att.clim,att.sim),
      #          max(att.sim.multi[[repName]][[tarName.1]],att.clim,att.sim))
      # plot(annAR1coeffList,att.sim.multi,ylim=ylim)
      # abline(h=att.clim,col='blue')
      # abline(h=att.sim,col='red')
      # abline(h=att.target,col='green')
      
      # browser()
      
      print(paste0(att.target,' ',att.sim.multi[[repName]][[tarName.1]][i]))
      
      annAR1coeff = annAR1coeffList[i]
      pp.annShuffle.out = pp.annShuffle(var.sim,times,annAR1coeff,seed,return.indices = T)
      #var.sim.targ = var.sim[pp.annShuffle.out$indices]

      var.sim = sim[[repName]][[tarName.1]][[perturb.varname]]$sim
      if (!is.matrix(var.sim)){var.sim=as.matrix(var.sim,ncol=1)}
      sim.new[[repName]][[tarName]][[perturb.varname]] = list(sim=var.sim[pp.annShuffle.out$indices,])
      for (other.varname in other.varnames){
        other.sim = sim[[repName]][[tarName.1]][[other.varname]]$sim
        if (!is.matrix(other.sim)){other.sim=as.matrix(other.sim,ncol=1)}
        other.sim.targ = other.sim[pp.annShuffle.out$indices,]
        sim.new[[repName]][[tarName]][[other.varname]] = list(sim=other.sim.targ)  
      }
      
      # clim_sim_shuffle = list(times=times)
      # clim_sim_shuffle[[perturb.varname]]=sim.new[[repName]][[tarName]][[perturb.varname]]$sim
      # print(calculateAttributes(clim_sim_shuffle,attPerturb))
      
      # PET.targ = PET.sim[o$indices]
      # 
      # PET.sim = sim[[repName]][[tarName.1]][['PET']]$sim
      
#      sim.new[[repName]][[tarName]][[var.2]] = sim.2[[repName.2]][[tarName.2]][[var.2]]
    }
  }

  sim.new$expSpace = expSpace.new
  sim.new$simDates = sim$simDates
  
  return(sim.new)
  
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
  #                                   P.sim=sim_stoch$Rep1$Target1$P$sim,
  #                                   target.dwellTime=dwellTime.target)
  # 
  
  
}

######################################################

