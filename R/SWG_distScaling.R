#################################

#' @include default_parameters.R

modelInfoList[["P-ann-distScaling"]] = list(simVar="P",
                                            timeStep='any',
                                       simPriority=1,
                                       npars=2,
                                       parNam=c("scale","shape"),
                                       minBound=c(1e-2,1e-2),
                                       maxBound=c(1e2,1e2))

# #################################

parManager.distScaling = function(parS, SWGparameterization, datInd, auxInfo=NULL){

  if (SWGparameterization=='ann'){
    scale <- parS['scale']
    shape <- parS['shape']
  }

  parTS = list(scale=scale,shape=shape)

  return(parTS)

}

#################################

SWGsim.distScaling = function(SWGpar,
                            nTimes=NULL,
                            randomTerm=NULL,
                            auxInfo){

  P = auxInfo$obs$P

  i.ww = which(P>0)
  dat = P[i.ww]

  pars.obs=gammaParsMLE2(dat=dat,wetThresh=0.00)

  CDFvals = pgamma(dat,shape=pars.obs$shape,scale=pars.obs$scale)

  P.new = P

  P.new[i.ww] = qgamma(CDFvals,shape=SWGpar$shape,scale=SWGpar$scale)

  return(P.new)

}

#################################
