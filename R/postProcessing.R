ppInfoList = list()

ppInfoList[['annVar']] = list(npars=1,
                              parNames=c('annVarFac'),
                              minBound=c(0.1),
                              maxBound=c(5))

ppInfoList[['scaleExtremes']] = list(scaleExtremesProb=0.99)

ppInfoList[['annCor']] = list(npars=1,
                              parNames=c('annAR1coeff'),
                              minBound=c(-0.5),
                              maxBound=c(0.99))

#################################

pp.annVar = function(sim,annVarFac,datInd){
  
  Pann = rep(NA,datInd$nyr)
  Pdaily = sim
  for(iy in 1:datInd$nyr){
    ind=datInd$i.yy[[iy]]
    Pann[iy] = sum(Pdaily[ind])
  } 
  
  meanPann = mean(Pann)
  Pann_new = meanPann + parS['annVarFac']*(Pann-meanPann)
  Pann_fac = Pann_new/Pann
  
  multSim=rep(NA,datInd$nTimes)
  for(iy in 1:datInd$nyr){
    ind=datInd$i.yy[[iy]]
    multSim[ind]=Pann_fac[iy]
  }
  multSim<-pmax(multSim,0)
  sim = Pdaily*multSim
  
  return(sim)
  
}

#################################

pp.scaleExtremes = function(sim,obs,prob){
  
  if (length(sim)!=length(obs)){browser()}
  
  fac = quantile(sim,prob) / quantile(obs,prob)
  nTop = floor(length(sim)*(1-prob))
  sortSim = sort(sim,decreasing = T,index.return=T)
  i = sortSim$ix[1:nTop]
  sortObs = sort(obs,decreasing = T)[1:nTop]
  
  sim[i] = fac*sortObs
  
  return(sim)
  
}

#################################

runPP = function(sim,obs,PPname,parS,datInd,randomTerm=NULL){
  
  if (PPname == 'annVar'){
    sim = pp.annVar(sim=sim,annVarFac=parS['annVarFac'],datInd=datInd)
  } else if (PPname == 'scaleExtremes'){
    sim = pp.scaleExtremes(sim=sim,obs=obs,prob=ppInfoList[[PPname]]$scaleExtremesProb)
  } else if (PPname == 'annCor'){
    sim = pp.annShuffle(P=sim,datInd=datInd,
                        annAR1coeff=ppInfoList[[PPname]]$annAR1coeff,
                        seed=randomTerm$seed)
  }
  
  return(sim)
  
}

###################

calcAC = function(x){
  nx = length(x)
  AC = cor(x[1:(nx-1)],x[2:nx])
  return(AC)
}

###################

shuffle = function(ar1coeff,sort.P.ann,index.P.ann,seed=1){
  
  set.seed(seed)
  ar1 = as.numeric(arima.sim(n=length(sort.P.ann),list(ar=ar1coeff)))
  
  rankAR1 = rank(ar1)
  
  P.ann.new = sort.P.ann[rankAR1]
  year.new=index.P.ann.long[rankAR1]
  
  P.ann.cor = calcAC(P.new) 
  
  out = list(P.ann.cor=P.ann.cor,P.new=P.new,year.new=year.new)

  return(out)
  
}

###################

shuffle_P = function(P,years,annAR1coeff,seed=1){
  
  browser()
  
  P.ann = c()
  for (y in 1:length(years)){
    year = years[y]
    keep = which(sim$simDates$year==year)
    P.ann[y] = sum(P[keep])
  }
  
  #  Nyears_long = length(years.long)
  Nyears_long = length(P.ann.long)
  
  s = sort(P.ann,index.return=T)
  sort.P.ann = s$x
  index.P.ann = s$ix
  
  o = shuffle(x=corVal,
              sort.P.ann=sort.P.ann,
              index.P.ann=index.P.ann.long,seed=seed,return.all = T)
  
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
  
  out[[label]]$P.ann.cor = o$P.ann.cor
  
  #################################
  
  P.new.daily = c()
  keepList = c()
  for (y in 1:Nyears_long){
    if (y%%1000==0){print(y)}
    year = years.long[1] + (year.new[y]-1)
    if (is.null(iyy)){
      keep = which(years.long==year)
    } else {
      keep = iyy[[year+1]]
    }
    keepList = c(keepList,keep)
  }
  
  P.new.daily = P.long[keepList]
  
  #      P.daily.changeCor[[label]] = P.new.daily
  
  out[[label]]$P.daily.changeCor = P.new.daily

  return(out)

}
#################################


