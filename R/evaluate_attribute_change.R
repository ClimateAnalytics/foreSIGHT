######################

systemModel_calcAtts <- function(data,          # data.frame with columns: year, month, day, *var1*, *var2* etc.
                                 systemArgs,    # list containing the arguments of simulateSystem
                                 metrics) {     # names of performance metrics (with units of the metrics)
  
  if (!is.null(systemArgs$vSel)){
    if (systemArgs$cSel=='mean'){
      data[[systemArgs$vSel]] = apply(data[[systemArgs$vSel]],1,mean)
    } else {
      data[[systemArgs$vSel]] = data[[systemArgs$vSel]][,systemArgs$cSel]
    }
  }
  
  systemPerformanceSim = calculateAttributes(climateData=data,attSel=metrics)
  
  systemPerformance = (systemPerformanceSim/systemArgs$attBase-1)*100
  
  return(systemPerformance)
}

#################################################################

#' @export
calcPerformanceAttributes = function(clim,sim,attSel,vSel=NULL,cSel=NULL){

  if (!is.null(cSel)){
    if (!is.null(vSel)){
      if (cSel=='mean'){
        clim[[vSel]] = apply(clim[[vSel]],1,mean)
      } else {
        clim[[vSel]] = clim[[vSel]][,cSel]
      }      
    } else {
      print('must enter vSel')
      return()
    }

  }

  attBase = calculateAttributes(clim,attSel)

  systemPerf <- runSystemModel(sim = sim,                     # simulation; the perturbed time series
                               systemModel = systemModel_calcAtts,      # the system model function
                               systemArgs = list(attBase=attBase,
                                                 vSel=vSel,
                                                 cSel=cSel),        # argument to the system model function
                               metrics = attSel)              # selected performance metrics

  return(systemPerf)
  
}

##################################################

#' @export
plotPerformanceAttributesOAT = function(clim,sim,attPerturb,attEval,Perf=NULL,vSel=NULL,cSel=NULL,
                                        ylim=NULL,
                                        cex.main=0.8,cex.xaxis=0.5,cex.yaxis=0.5){

  if (is.null(Perf)){
    Perf = calcPerformanceAttributes(clim=clim,sim=sim,attSel=attEval,vSel=vSel,cSel=cSel)
  }

  for (att in names(Perf)){
    o = plotPerformanceOAT(Perf, sim, metric=att,attSel=attPerturb,returnPlotData=T)
    plotPerformanceOAT.baseR(plotData=o$plotData,sim=sim,metric=att,targetVal = o$targetVal,
                             attSel=attPerturb,
                             ylim=ylim,
                             cex.main=cex.main,cex.xaxis=cex.xaxis,cex.yaxis=cex.yaxis)
  } 

}

##################################################

plotPerformanceOAT.baseR = function(plotData,sim,metric,attSel,
                                    targetVal,
                                    ylim=NULL,baseSettings=list(),
                                    cex.main=0.8,cex.xaxis=0.5,cex.yaxis=0.5){
  
  if(is.null(baseSettings$bias_base_thresh)){baseSettings$bias_base_thresh = 20}
  if(is.null(baseSettings$slope_thresh)){baseSettings$slope_thresh = 1.5}
  
  colMed = 'black'; lwdMed = 1
  colShade='lightgrey'
  colHold = 'green'; ltyHold = 2; lwdHold = 3
  colTied = 'cyan'; ltyTied = 2; lwdTied = 3
  colPert = 'blue'; ltyPert = 2; lwdPert = 3
  colZero = 'black'; ltyZero = 3; lwdZero=0.5 
  lwd=1
  
  # # determine target values associated with OAT perturbations  
  # targetMat = sim$expSpace$targetMat
  # if (metric %in% colnames(targetMat)){
  #   targetVal =targetMat[iInd,metric]
  #   targetVal = (targetVal-1)*100
  # } else {
  #   targetVal = NULL
  # }
  
  # determine median and upper and lower limits
  m = 1
  x = plotData[[m]][,1]
  x = (x-1)*100
  med = plotData[[m]][,2]
  if (dim(plotData[[m]])[2]==5){
    lo = plotData[[m]][,3]
    hi = plotData[[m]][,4]
  } else {
    lo=NULL
    hi=NULL
  }
  
  if (is.null(ylim)){
    yMin = min(med,lo,hi,-10)
    yMax = max(med,lo,hi,10)
    ylim = c(yMin,yMax)
  }
  
  plot(x=x,y=med,type='o',ylim=ylim,xaxs='i',xlab='',ylab='',col=colMed)
  if (!is.null(lo)){
    polygon(c(rev(x), x), c(rev(hi), lo), col = colShade, border = NA)
  }
  lines(x,med,type='l',col=colMed,lwd=lwdMed)
  box()
  if(!is.null(targetVal)){
    if (metric==attSel){
      col=colPert
      lty=ltyPert
      lwd=lwdPert
    } else {
      col=colHold
      lty=ltyHold
      lwd=lwdHold          
    }
    lines(x,targetVal,col=col,lty=lty,lwd=lwd)
  }
  points(x,med,col=colMed,lwd=lwdMed)
  
  abline(h=0,lty=ltyZero,lwd=lwdZero,col=colZero)
  abline(v=0,lty=ltyZero,lwd=lwdZero,col=colZero)
  
  bias_base = med[x==0]
  bias_base_hi = abs(bias_base) > baseSettings$bias_base_thresh
  if (bias_base_hi){points(x=0,bias_base,col='red',pch=4,cex=2,lwd=2)}
  
  mod = lm(med~x)
  slope = mod$coefficients[2] 
  inflated_response = abs(slope)>baseSettings$slope_thresh  
  if (inflated_response){lines(x,med,col='red',lwd=2)}
  
  title_str = metric
  if (bias_base_hi){title_str=paste0(title_str,' B')}
  if (inflated_response){title_str=paste0(title_str,' I')}
  title(title_str,cex.main=cex.main)
  
  #attribute = unique(plotData[[m]][,'attribute'])
  #mtext(side=1,text=attribute,line = 2,cex = 0.7)
  
  # mtext(side=1,text=paste0('D ',attPerturb,' (%)'),line = 2,cex = cex.xaxis)
  # mtext(side=2,text=paste0('D ',metric,' (%)'),line = 2,cex = cex.yaxis)
  
  mtext(side=1,text='Change pert att (%)',line = 2,cex = cex.xaxis)
  mtext(side=2,text='Change att (%)',line = 2,cex = cex.yaxis)
  
}


