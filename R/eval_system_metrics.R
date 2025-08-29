#' @export
evaluate_system_metrics = function(sim,clim,systemModel,systemArgs,metrics,obs_metrics=NULL,varNames=NULL){
  
  expSpace = sim$expSpace
  
  #########
  # determine target number that represent baseline unperturbed climate 
  attSel = colnames(sim$expSpace$targetMat)
  varType <- vapply(attSel, FUN = get.attribute.varType, FUN.VALUE = character(1), USE.NAMES = FALSE)
  targetType <- vapply(varType, FUN = get.target.type, FUN.VALUE = character(1), USE.NAMES = FALSE)
  
  baseVal = rep(NA,length(attSel))
  baseVal[targetType=='diff'] = 0
  baseVal[targetType=='frac'] = 1
  
  targets = expSpace$targetMat
  b = which(apply(targets==baseVal,1,FUN =all))
  
  if (length(b)==0){
    print('no baseline scenario')
    return()
  }
  
  #########
  # strip other targets from sim
  simBase = sim
  numReps = length(which(grepl('Rep',names(sim))))
  numTars = length(names(sim[[1]]))
  for (r in 1:numReps){
    repName = paste0('Rep',r)
    simBase[[repName]] = NULL
    simBase[[repName]][['Target1']] = sim[[repName]][[paste0('Target',b)]] 
  }
  
  #########
  # performance for baseline unperturbed climate  
  systemPerf_base <- runSystemModel(sim = simBase,                     # simulation; the perturbed time series
                                    systemModel = systemModel,      # the system model function
                                    systemArgs = systemArgs,        # argument to the system model function
                                    metrics = metrics,
                                    varNames=varNames)              # selected performance metrics 
  
  # performance using observed climate 
  systemPerf_obsClim = systemModel(data = clim, systemArgs = systemArgs, metrics = metrics)
  
  #############
  
  return(list(systemPerf_base=systemPerf_base,
              systemPerf_obsClim=systemPerf_obsClim))
  
}

#' Draws a boxplot with the whiskers at specified probability limits 
#'
#' \code{boxplot_prob} Draws a boxplot with the whiskers at the probability limits, provided by \code{whiskersProb}.
#' @param xin a vector, matrix or dataframe; data to be plotted
#' @param whiskersProb a vector of length 2; min and max probability limits
#' @param at a vector; specifying x coordinates for boxes
#' @return The function returns a boxplot. 
#' @export
boxplot_prob =function(xin,whiskersProb=c(0.025,0.975),at=NULL,...){
  # Draws a boxplot with the whiskers at the probability limits, provided by whiskersProb
  # x can be vector or data.frame
  # whiskers Prob values are probabilities for the lower and upper whisker (respectively)
  
  x = data.frame(xin)
  x.stats=graphics::boxplot(x,plot=FALSE)
  
  if ( is.data.frame(x) || is.matrix(x)) {
    x.stats$out = c()
    x.stats$group = c()
    y=x
    for (j in 1:ncol(x)) {
      y=sort(x[,j])
      if (length(stats::na.omit(y))!=0) {
        x.stats$stats[1,j]=stats::quantile(y,prob=whiskersProb[1])
        x.stats$stats[5,j]=stats::quantile(y,prob=whiskersProb[2])
      }
      for (i in 1:length(stats::na.omit(y))){
        if (y[i]<x.stats$stats[1,j]|y[i]>x.stats$stats[5,j]){
          
          x.stats$out = c(x.stats$out,y[i])
          x.stats$group = c(x.stats$group,j)
          
        }
      }
      
    }
  } else if (is.vector(x)) { # Assume x is a vector
    x.stats$out = c()
    x.stats$group = c()
    y=sort(x)
    if (length(stats::na.omit(y))!=0) {
      x.stats$stats[1,1]=stats::quantile(y,prob=whiskersProb[1])
      x.stats$stats[5,1]=stats::quantile(y,prob=whiskersProb[2])
    }
    for (i in 1:length(stats::na.omit(y))){
      if (y[i]<x.stats$stats[1]|y[i]>x.stats$stats[5]){
        
        x.stats$out = c(x.stats$out,y[i])
        x.stats$group = c(x.stats$group,1)
        
      }
    }
  } else
  {
    print("type of x is not supported in boxplot.ext")
    return()
  }
  
  graphics::bxp(z=x.stats,at=at,...)
  
}

