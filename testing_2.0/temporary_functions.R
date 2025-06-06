mvFunc_cor = function(data.1,data.2){
  return(cor(data.1,data.2,use='pairwise.complete.obs'))
}

mvFunc_avgWetDay = function(data.1,data.2){
  return(mean(data.1[data.2>0],na.rm=T))
}

mvFunc_sdWetDay = function(data.1,data.2){
  return(sd(data.1[data.2>0],na.rm=T))
}

mvFunc_avgDryDay = function(data.1,data.2){
  return(mean(data.1[data.2==0],na.rm=T))
}

mvFunc_sdDryDay = function(data.1,data.2){
  return(sd(data.1[data.2==0],na.rm=T))
}

func_sd = function(data){
  return(sd(data,na.rm=T))
}

func_xP90 = function(data){
  P90 = quantile(data,probs = 0.9,na.rm=T,names=F)
  if(is.na(P90)){browser()}
  return(P90)
}

mvFunc_xP90WetDay = function(data.1,data.2){
  return(quantile(data.1[data.2>0],probs=0.9,na.rm=T,names=F))
}

mvFunc_xP90DryDay = function(data.1,data.2){
  return(quantile(data.1[data.2==0],probs=0.9,na.rm=T,names=F))
}


