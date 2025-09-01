######################

systemModel_calcAtts <- function(data,          # list with climate data and times
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
# calculate changes in attributes over all targets and replicates 
calcPerformanceAttributes = function(clim,      # reference climate
                                     sim,       # simulated climates from generateScenarios()
                                     attSel,    # attributes to calculate
                                     vSel=NULL, # variable name for selection site/aggregation
                                     cSel=NULL)  # how to summarize multisite data - select site number or 'mean' for average 
  {

  # convert multiste data to single time series
  if (!is.null(cSel)){
    if (!is.null(vSel)){
      if (cSel=='mean'){
        clim[[vSel]] = apply(clim[[vSel]],1,mean)
      } else {
        clim[[vSel]] = clim[[vSel]][,cSel]
      }      
    } else {
      stop('must enter vSel')
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

#' Plots changes in attributes for a specified perturbed attribute
#'
#' \code{plotPerformanceAttributesOAT} plots OAT changes in attributes based on simulated climate object, for a given perturbed attribute.
#' 
#' @param clim list; reference climate data. \cr 
#' @param sim list; perturbed climate from \code{generateScenarios()} \cr 
#' @param attPerturb string; name of perturbed attribute \cr
#' @param attEval string; name of attribute that will be evaluated \cr
#' @param vSel string; variable name for selection site/aggregation \cr
#' @param cSel integer or string 'mean'; how to summarize multisite data - select site number or 'mean' for average \cr 
#' @param ylim numeric vector of length 2; min and max y limits for plotting \cr
#' @param cex.main number; size for title \cr
#' @param cex.xaxis number; size for x-axis \cr
#' @param cex.yaxis number; size for y-axis \cr
#' @return The function returns a single plot showing changes in attribute \code{attEval} for changes in perturbed attribute \code{attPerturb} \cr  
#' @examples 
#' # XXXXXXXXXX
#' @export
plotPerformanceAttributesOAT = function(clim,
                                        sim,
                                        attPerturb,
                                        attEval,
#                                        Perf=NULL,
                                        vSel=NULL,cSel=NULL,
                                        ylim=NULL,
                                        cex.main=0.8,cex.xaxis=0.5,cex.yaxis=0.5){

#  if (is.null(Perf)){
   Perf = calcPerformanceAttributes(clim=clim,sim=sim,attSel=attEval,vSel=vSel,cSel=cSel)
#  }

  for (att in names(Perf)){
    o = plotPerformanceOAT(Perf, sim, metric=att,attSel=attPerturb,returnPlotData=T)
    plotPerformanceOAT.baseR(plotData=o$plotData,sim=sim,metric=att,targetVal = o$targetVal,
                             attSel=attPerturb,
                             ylim=ylim,
                             cex.main=cex.main,cex.xaxis=cex.xaxis,cex.yaxis=cex.yaxis)
  } 

}

##################################################
# function called from plotPerformanceAttributesOAT()
# code for plotting changes in attributes for a given perturbed attribute
plotPerformanceOAT.baseR = function(plotData,      # list containing changes in attributes (calculated from plotPerformanceOAT)
                                    sim,           # sim object from generateScenarios()
                                    metric,        # attribute name
                                    attSel,        # perturbed attribute
                                    targetVal,     # target value of attribute
                                    ylim=NULL,     # y limit
                                    baseSettings=list(),  # list containing 'bias_base_thresh' for threshold for bias in baseline performance (%), and 'slope_thresh' for threshold for slope of relationship between perturbed and plotted attribute
                                    cex.main=0.8,cex.xaxis=0.5,cex.yaxis=0.5) # font size for title and axis labels
  {
  
  if(is.null(baseSettings$bias_base_thresh)){baseSettings$bias_base_thresh = 20}
  if(is.null(baseSettings$slope_thresh)){baseSettings$slope_thresh = 1.5}
  
  colMed = 'black'; lwdMed = 1
  colShade='lightgrey'
  colHold = 'green'; ltyHold = 2; lwdHold = 3
  colTied = 'cyan'; ltyTied = 2; lwdTied = 3
  colPert = 'blue'; ltyPert = 2; lwdPert = 3
  colZero = 'black'; ltyZero = 3; lwdZero=0.5 
  colPoorPerformance = 'red'
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
    graphics::polygon(c(rev(x), x), c(rev(hi), lo), col = colShade, border = NA)
  }
  graphics::lines(x,med,type='l',col=colMed,lwd=lwdMed)
  graphics::box()
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
    graphics::lines(x,targetVal,col=col,lty=lty,lwd=lwd)
  }
  graphics::points(x,med,col=colMed,lwd=lwdMed)
  
  graphics::abline(h=0,lty=ltyZero,lwd=lwdZero,col=colZero)
  graphics::abline(v=0,lty=ltyZero,lwd=lwdZero,col=colZero)
  
  bias_base = med[x==0]
  bias_base_hi = abs(bias_base) > baseSettings$bias_base_thresh
  if (bias_base_hi){graphics::points(x=0,bias_base,col=colPoorPerformance,pch=4,cex=2,lwd=2)}
  
  mod = stats::lm(med~x)
  slope = mod$coefficients[2] 
  inflated_response = abs(slope)>baseSettings$slope_thresh  
  if (inflated_response){graphics::lines(x,med,col=colPoorPerformance,lwd=2)}
  
  title_str = metric
  if (bias_base_hi){title_str=paste0(title_str,' B')}
  if (inflated_response){title_str=paste0(title_str,' I')}
  graphics::title(title_str,cex.main=cex.main)
  
  #attribute = unique(plotData[[m]][,'attribute'])
  #mtext(side=1,text=attribute,line = 2,cex = 0.7)
  
  # mtext(side=1,text=paste0('D ',attPerturb,' (%)'),line = 2,cex = cex.xaxis)
  # mtext(side=2,text=paste0('D ',metric,' (%)'),line = 2,cex = cex.yaxis)
  
  graphics::mtext(side=1,text='Change pert att (%)',line = 2,cex = cex.xaxis)
  graphics::mtext(side=2,text='Change att (%)',line = 2,cex = cex.yaxis)
  
}


