ppInfoList = list()

ppInfoList[['annVar']] = list(npars=1,
                              parNam=c('annVarFac'),
                              minBound=c(0.1),
                             maxBound=c(5))

ppInfoList[['scaleExtremesAll']] = list(npars=0,
                                        scaleExtremesProb=0.99)

ppInfoList[['scaleExtremesSeas']] = list(npars=0,
                                         scaleExtremesProb=0.99)

ppInfoList[['annCor']] = list(npars=1,
                              parNam=c('annAR1coeff'),
                              minBound=c(-0.5),
                              maxBound=c(0.99))

#################################

runPP = function(sim,obs,PPname,parS,datInd,randomTerm=NULL){
  
  if (PPname == 'annVar'){
    sim = pp.annVar(sim=sim,annVarFac=parS['annVarFac'],datInd=datInd)
  } else if (PPname == 'scaleExtremesAll'){
    sim = pp.scaleExtremes(sim=sim,obs=obs,prob=ppInfoList[[PPname]]$scaleExtremesProb,strat='all',datInd=datInd)
  } else if (PPname == 'scaleExtremesSeas'){
    sim = pp.scaleExtremes(sim=sim,obs=obs,prob=ppInfoList[[PPname]]$scaleExtremesProb,strat='seas',datInd=datInd)
  } else if (PPname == 'annCor'){
    sim = pp.annShuffle(P=sim,datInd=datInd,
                        annAR1coeff=parS['annAR1coeff'],
                        seed=randomTerm$seed)
  }
  
  return(sim)
  
}

#################################

pp.annVar = function(sim,annVarFac,datInd){
  
  Pann = rep(NA,datInd$nyr)
  Pdaily = sim
  for(iy in 1:datInd$nyr){
    ind=datInd$i.yy[[iy]]
    Pann[iy] = sum(Pdaily[ind])
  } 

Pann[Pann==0] = 0.00001

  meanPann = mean(Pann)
  Pann_new = meanPann + annVarFac*(Pann-meanPann)
  Pann_fac = Pann_new/Pann
  
  multSim=rep(NA,datInd$nTimes)
  for(iy in 1:datInd$nyr){
    ind=datInd$i.yy[[iy]]
    multSim[ind]=Pann_fac[iy]
  }
  multSim<-pmax(multSim,0)
  sim = Pdaily*multSim
 
   if(any(is.na(sim))){
  browser()
}
 
  return(sim)
  
}

#################################

pp.scaleExtremes = function(sim,obs,prob,strat,datInd){
  
  if (length(sim)!=length(obs)){browser()}
  
  # fac = quantile(sim,prob) / quantile(obs,prob)
  # nTop = floor(length(sim)*(1-prob))
  # sortSim = sort(sim,decreasing = T,index.return=T)
  # i = sortSim$ix[1:nTop]
  # sortObs = sort(obs,decreasing = T)[1:nTop]
  # sim[i] = fac*sortObs

  if (strat=='all'){
    N = 1
  } else if (strat=='seas'){
    N = 4
  }
  
  for (s in 1:N){
    if (strat=='all'){
      keep = 1:datInd$nTimes
    } else if (strat=='seas'){
      keep = datInd$i.ss[[s]]
    }
  if (any(is.na(sim))){print('sim has na')}
  if (any(is.na(sim[keep]))){print('sim[keep] has na')}
    fac = quantile(sim[keep],prob) / quantile(obs[keep],prob)
    nTop = floor(length(sim[keep])*(1-prob))
    sortSim = sort(sim[keep],decreasing = T,index.return=T)
    i = keep[sortSim$ix[1:nTop]]
    sortObs = sort(obs[keep],decreasing = T)[1:nTop]
    sim[i] = fac*sortObs
  }
  
  return(sim)
  
}


#   # for (s in 1:4){
#   #   keep = modelEnv$P_modelEnv$datInd$i.ss[[s]]
#   #   fac = quantile(sim$sim[keep],0.99) / quantile(obs$P[keep],0.99)
#   #   nTop = floor(length(sim$sim[keep])/100)
#   #   tmp.sortSim = sort(sim$sim[keep],decreasing = T,index.return=T)
#   #   i = keep[tmp.sortSim$ix[1:nTop]]
#   #   sortObs = sort(obs$P[keep],decreasing = T)[1:nTop]
#   #   sim$sim[i] = fac*sortObs
#   # }
###################

calcAC = function(x){
  nx = length(x)
  AC = cor(x[1:(nx-1)],x[2:nx])
  return(AC)
}

###################

shuffle = function(annAR1coeff,sort.P.ann,index.P.ann,seed=1){

  set.seed(seed)
  ar1 = suppressWarnings(as.numeric(arima.sim(n=length(sort.P.ann),list(ar=annAR1coeff))))

  rankAR1 = rank(ar1)

  P.ann.new = sort.P.ann[rankAR1]
  year.new=index.P.ann[rankAR1]

  P.ann.cor = calcAC(P.ann.new)

  out = list(P.ann.cor=P.ann.cor,
             P.new=P.ann.new,
             year.new=year.new)

  return(out)

}

###################

pp.annShuffle = function(P,times,annAR1coeff,seed=1,iyy=NULL,return.indices=F){

  # print(annAR1coeff)
  
  years.all = as.integer(format(times,'%Y'))
  years = unique(years.all)
  
  P.ann = c()
  for (y in 1:length(years)){
    year = years[y]
    keep = which(years.all==year)
    P.ann[y] = sum(P[keep])
  }
  
  s = sort(P.ann,index.return=T)
  sort.P.ann = s$x
  index.P.ann = s$ix

  o = shuffle(annAR1coeff=annAR1coeff,
              sort.P.ann=sort.P.ann,
              index.P.ann=index.P.ann,
              seed=seed)

  P.ann.cor=o$P.ann.cor
  P.new=o$P.new
  year.new=o$year.new

  # ar1 = as.numeric(arima.sim(n=Nyears_long,list(ar=corVal)))
  #
  # rankAR1 = rank(ar1)
  #
  # P.new = sort.P.ann.long[rankAR1]
  # year.new=index.P.ann.long[rankAR1]
  #
  # out[[label]]$P.ann.cor = calcAC(P.new)

  # out[[label]]$P.ann.cor = o$P.ann.cor

  #################################

  P.new.daily = c()
  keepList = c()
  for (y in 1:length(years)){
    year = years[1] + (year.new[y]-1)
    if (is.null(iyy)){
      keep = which(years.all==year)
    } else {
      keep = iyy[[year+1]]
    }
    keepList = c(keepList,keep)
  }

  P.new.daily = P[keepList]

  #      P.daily.changeCor[[label]] = P.new.daily

  # climSim = list(times=datInd$times,P=P.new.daily)
  # print(calculateAttributes(climSim,'P_day_all_tot_cor'))
  
  if (!return.indices){
    return(P.new.daily)
  } else {
    return(list(P.new=P.new.daily,
                indices=keepList))
  }

}
#################################


