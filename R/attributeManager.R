#######################################
##      ATTRIBUTE MANAGER            ##
#######################################

#CONTAINS
  #functions of form "func_XXX", used to calculate attributes
  #attribute.calculator() - calculate values of attributes
  #attribute.calculator.setup() - calculate arguments used in attribute.calculator() based on attribute names
  #attribute.info.check() - get targetType, varType and identify any invalid model selections
    #check.attribute.model.combo() - check if any attribute-model combos are invalid
    #get.attribute.info() - Identify invalid models
    #get.target.type() - treated as fractions, percent or abs value
    #get.attribute.varType() - "P", "Temp" OR ...


#---------------------------------------------------------------------------------------

####################################
# Functions used to calculate attributes.
# Each has input of "data" and optional "attArgs"
# Note that custom functions starting with "func_" can also be specified

#' Calculates total of time series
#' @param data is a vector, representing a time series
#' @export
func_tot = function(data) sum(data)

#' Calculates seasonality ratio
#' @param data is a vector, representing a time series
#' @param attArgs is a list, with attArgs$indexWet corresponding to wet season and attArgs$indexDry dry season
# seasonality ratio
#' @export
func_seasRatio = function(data,attArgs){
  Pdry = sum(data=data[attArgs$indexDry])
  Pwet = sum(data=data[attArgs$indexWet])
  Pdry = max(Pdry,0.001)
  seasRatio = Pwet/Pdry
  return(seasRatio)
}

#' Calculates number of wet days (above threshold)
#' @param data is a vector, representing a time series
#' @param attArgs is a list, with attArgs$threshold denoting the threshold
#' @export
func_nWet = function(data,attArgs) get.nwet(data=data,threshold=attArgs$threshold)

#' Calculates maximum dry spell duration (below threshold)
#' @param data is a vector, representing a time series
#' @param attArgs is a list, with attArgs$threshold denoting the threshold
#' @export
func_maxDSD = function(data,attArgs) get.spell.lengths.max(data=data,thresh=attArgs$threshold,type="dry")

#' Calculates maximum wet spell duration (above threshold)
#' @param data is a vector, representing a time series
#' @param attArgs is a list, with attArgs$threshold denoting the threshold
#' @export
func_maxWSD = function(data,attArgs) get.spell.lengths.max(data=data,thresh=attArgs$threshold,type="wet")

#' Calculates average dry spell duration (below threshold)
#' @param data is a vector, representing a time series
#' @param attArgs is a list, with attArgs$threshold denoting the threshold
#' @export
func_avgDSD = function(data,attArgs) mean(get.spell.lengths(data=data,thresh=attArgs$threshold,type="dry"),na.rm=TRUE)

#' Calculates average wet spell duration (below threshold)
#' @param data is a vector, representing a time series
#' @param attArgs is a list, with attArgs$threshold denoting the threshold
#' @export
func_avgWSD = function(data,attArgs) mean(get.spell.lengths(data=data,thresh=attArgs$threshold,type="wet"),na.rm=TRUE)

#' Calculates average rainfall on wet days (above threshold)
#' @param data is a vector, representing a time series
#' @param attArgs is a list, with attArgs$threshold denoting the threshold
#' @export
func_dyWet = function(data,attArgs) get.wet.average(data=data,threshold=attArgs$threshold)

#' Calculates the number of days above a threshold (often used for temperature)
#' @param data is a vector, representing a time series
#' @param attArgs is a list, with attArgs$threshold denoting the threshold
#' @export
func_R = function(data,attArgs) get.nwet(data=data,threshold=attArgs$threshold)

#' Calculates a quantile value
#' @param data is a vector, representing a time series
#' @param attArgs is a list, with attArgs$quant denoting the probability of the quantile
#' @export
func_P = function(data,attArgs) get.quantile(data=data,quant=attArgs$quant)

#' Calculates average of time series
#' @param data is a vector, representing a time series
#' @export
func_avg = function(data) mean(data,na.rm=T)

# inter-quantile range
#' Calculates the inter-quantile range
#' @param data is a vector, representing a time series
#' @param attArgs is a list, with attArgs$lim denoting the probability limit width
#' @export
func_rng = function(data,attArgs) get.quantile.rng(data=data,lim=attArgs$lim)

#' Calculates the growing season length
#' @param data is a vector, representing a time series
#' @export
func_GSL = function(data) GSLcalc(x=data)

#' Calculates the cold season length
#' @param data is a vector, representing a time series
#' @export
func_CSL = function(data) CSLcalc(x=data)

#' Calculates the number of frost days
#' @param data is a vector, representing a time series
#' @export
func_F0 = function(data) F0calc(x=data) # could be made generic

#' Calculates the day of year corresponding to the wettest 6 months
#' @param data is a vector, representing a time series
#' @param attArgs is a list, with attArgs$doy denoting the day of year for each value in the time series
#' @export
func_wettest6monPeakDay = function(data,attArgs=NULL){
#  if (is.null(attArgs$seas)){
#    seas = calc_meanClimDaily_dayOfYearWindow(obs=data,doy=attArgs$doy,inc=91)
#  } else {
    seas = attArgs$seas
#  }
  i = stats::median(which(seas==max(seas)))
  #  print(i)
  #  return(i-180)
}

#' Calculates the ratio of wet season to dry season rainfall, based on wettest6monPeakDay
#' @param data is a vector, representing a time series
#' @param attArgs is a list, with attArgs$doy denoting the day of year for each value in the time series
#' @export
func_wettest6monSeasRatio = function(data,attArgs=NULL){
#  if (is.null(attArgs$seas)){
#    seas = calc_meanClimDaily_dayOfYearWindow(obs=data,doy=attArgs$doy,inc=91)
#  } else {
    seas = attArgs$seas
#  }
  iwet = stats::median(which(seas==max(seas)))
  idry = stats::median(which(seas==min(seas)))
#  wettest6monSeasRatio = seas[iwet]/seas[idry]
  seas_iwet = seas[iwet]
  seas_idry = seas[idry]
  wettest6monSeasRatio = seas_iwet/seas_idry
  if ((seas_idry==0.)|(wettest6monSeasRatio>100.)){wettest6monSeasRatio=100.}
  return(wettest6monSeasRatio)
}

# for each doy calculate which dates have that doy, store results in matrix
calc_keepMat = function(doy){
  keepMat = matrix(nrow=365,ncol=length(doy)/365)
  for (n in 1:365){
    keepMat[n,] = which(doy==n)
  }
  return(keepMat)
}

calc_mean_day_clim = function(obs,keepMat){
  mean_day_clim = c()
  for (i in 1:365){
    mean_day_clim[i] = mean(obs[keepMat[i,]])
  }
  return(mean_day_clim)
}

# Calculates the seasonal pattern (i.e. climatological mean)
calc_meanClimDaily_dayOfYearWindow = function(obs,  # vector representing a time series
                                              dates=NULL,
                                              keepMat=NULL, # matrix indexing days of year
                                              inc) #the half-window size used in moving average
{
  if (is.null(keepMat)){
    doy = as.integer(format(dates,'%j'))
    keepMat = calc_keepMat(doy)
  }
  mean_day_clim = calc_mean_day_clim_cpp(obs,keepMat)
  indicesRM = c( (365-inc+1):365 , 1:365, 1:inc )
  run_mean_day_clim = ma(mean_day_clim[indicesRM],n=(2*inc+1))
  return(run_mean_day_clim[(inc+1):(inc+365)])
}

# Calculates the seasonal pattern (i.e. climatological mean) over all dates (not just single year)
calc_meanClimDaily_dayOfYearWindow_allDates = function(obs,  # vector representing a time series
                                              dates,
                                              inc) #the half-window size used in moving average
{

  obsClim = calc_meanClimDaily_dayOfYearWindow (obs=obs,dates = dates,inc = inc)
  obsClim[366] = obsClim[365]
  doy = as.integer(format(dates,'%j'))
  obsClimAll = obsClim[doy]

  return(obsClimAll)
}

####################################
#ATTRIBUTE CALCULATOR FUNCTION
attribute.calculator<-function(attSel=NULL,         #list of evaluated attribute names
                               data=NULL,           #timeseries data
                               datInd=NULL,         #dat indices and properties (e.g. datInd$nyr, datInd$i.yy)
                               attInfo=NULL         #optional saved list of attribute information (from attribute.calculator.setup)
){

  if (!is.null(attInfo$attCalcInfo)){
    attCalcInfo = attInfo$attCalcInfo
  } else {
    attCalcInfo = attribute.calculator.setup(attSel,datInd)
  }

  if (any(c("P_ann_wettest6monSeasRatio","P_ann_wettest6monPeakDay")%in%attSel)){
    seas = calc_meanClimDaily_dayOfYearWindow(obs=data,keepMat=attCalcInfo[["P_ann_wettest6monSeasRatio"]]$attArgs$keepMat,inc=91)
    attCalcInfo[["P_ann_wettest6monSeasRatio"]]$attArgs$seas = attCalcInfo[["P_ann_wettest6monPeakDay"]]$attArgs$seas = seas
  }

  out = list()
  for (att in attSel){
    if (is.null(attCalcInfo[[att]]$opName)){ # case where there is no operator in attribute name
      if (is.null(dim(data))){
        out[[att]] = extractor(func=attCalcInfo[[att]]$func,
                               data=data,
                               indx=attCalcInfo[[att]]$indx,
                               attArgs=attCalcInfo[[att]]$attArgs)
      } else {
        out[[att]] = apply(X=data,MARGIN=2,FUN=extractor,
                           func=attCalcInfo[[att]]$func,
                           indx=attCalcInfo[[att]]$indx,
                           attArgs=attCalcInfo[[att]]$attArgs)
      }
    } else if (attCalcInfo[[att]]$opName%in%c('m','m10yrBlock')){ # mean of values calculated in each year
      if (is.null(dim(data))){
        out[[att]] = extractor.summaryMean(func=attCalcInfo[[att]]$func,
                                           data=data,
                                           indx=attCalcInfo[[att]]$indx,
                                           attArgs=attCalcInfo[[att]]$attArgs)
      } else {
        out[[att]] = apply(X=data,MARGIN=2,FUN=extractor.summaryMean,
                           func=attCalcInfo[[att]]$func,
                           indx=attCalcInfo[[att]]$indx,
                           attArgs=attCalcInfo[[att]]$attArgs)
      }
    } else if (attCalcInfo[[att]]$opName=='sd'){ # sd of values calculated in each year
      if (is.null(dim(data))){
        out[[att]] = extractor.summarySD(func=attCalcInfo[[att]]$func,
                                           data=data,
                                           indx=attCalcInfo[[att]]$indx,
                                           attArgs=attCalcInfo[[att]]$attArgs)
      } else {
        out[[att]] = apply(X=data,MARGIN=2,FUN=extractor.summarySD,
                           func=attCalcInfo[[att]]$func,
                           indx=attCalcInfo[[att]]$indx,
                           attArgs=attCalcInfo[[att]]$attArgs)
      }
####### NOTE: calling the following separately is inefficient (annual totals calculated for each)
    } else if (attCalcInfo[[att]]$opName=='cor'){ # correlation between years
      if (is.null(dim(data))){
        out[[att]] = extractor.summaryCor(func=attCalcInfo[[att]]$func,
                                         data=data,
                                         indx=attCalcInfo[[att]]$indx,
                                         attArgs=attCalcInfo[[att]]$attArgs)
      } else {
        out[[att]] = apply(X=data,MARGIN=2,FUN=extractor.summaryCor,
                           func=attCalcInfo[[att]]$func,
                           indx=attCalcInfo[[att]]$indx,
                           attArgs=attCalcInfo[[att]]$attArgs)
      }
    } else if (attCalcInfo[[att]]$opName=='corSOI'){ # 
      if (is.null(dim(data))){
        out[[att]] = extractor.summaryCorSOI(func=attCalcInfo[[att]]$func,
                                          data=data,
                                          indx=attCalcInfo[[att]]$indx,
                                          attArgs=attCalcInfo[[att]]$attArgs)
      } else {
        out[[att]] = apply(X=data,MARGIN=2,FUN=extractor.summaryCorSOI,
                           func=attCalcInfo[[att]]$func,
                           indx=attCalcInfo[[att]]$indx,
                           attArgs=attCalcInfo[[att]]$attArgs)
      }
    } else if (attCalcInfo[[att]]$opName=='dwellTime'){ # sd of values calculated in each year
      if (is.null(dim(data))){
        out[[att]] = extractor.summaryDwellTime(func=attCalcInfo[[att]]$func,
                                          data=data,
                                          indx=attCalcInfo[[att]]$indx,
                                          attArgs=attCalcInfo[[att]]$attArgs)
      } else {
        out[[att]] = apply(X=data,MARGIN=2,FUN=extractor.summaryDwellTime,
                           func=attCalcInfo[[att]]$func,
                           indx=attCalcInfo[[att]]$indx,
                           attArgs=attCalcInfo[[att]]$attArgs)
      }
    } else if (attCalcInfo[[att]]$opName=='range90'){ # sd of values calculated in each year
      if (is.null(dim(data))){
        out[[att]] = extractor.summaryRange90(func=attCalcInfo[[att]]$func,
                                                data=data,
                                                indx=attCalcInfo[[att]]$indx,
                                                attArgs=attCalcInfo[[att]]$attArgs)
      } else {
        out[[att]] = apply(X=data,MARGIN=2,FUN=extractor.summaryRange90,
                           func=attCalcInfo[[att]]$func,
                           indx=attCalcInfo[[att]]$indx,
                           attArgs=attCalcInfo[[att]]$attArgs)
      }
    } else if (attCalcInfo[[att]]$opName%in%c('max3yr','max5yr')){ # maximum of 3 or 5 year values
      if (is.null(dim(data))){
        out[[att]] = extractor.summaryMax(func=attCalcInfo[[att]]$func,
                                          data=data,
                                          indx=attCalcInfo[[att]]$indx,
                                          attArgs=attCalcInfo[[att]]$attArgs)
      } else {
        out[[att]] = apply(X=data,MARGIN=2,FUN=extractor.summaryMax,
                           func=attCalcInfo[[att]]$func,
                           indx=attCalcInfo[[att]]$indx,
                           attArgs=attCalcInfo[[att]]$attArgs)
      }
    } else if (attCalcInfo[[att]]$opName%in%c('min3yr','min5yr')){ # maximum of 3 or 5 year values
      if (is.null(dim(data))){
        out[[att]] = extractor.summaryMin(func=attCalcInfo[[att]]$func,
                                          data=data,
                                          indx=attCalcInfo[[att]]$indx,
                                          attArgs=attCalcInfo[[att]]$attArgs)
      } else {
        out[[att]] = apply(X=data,MARGIN=2,FUN=extractor.summaryMin,
                           func=attCalcInfo[[att]]$func,
                           indx=attCalcInfo[[att]]$indx,
                           attArgs=attCalcInfo[[att]]$attArgs)
      }
    } else if (attCalcInfo[[att]]$opName%in%c('p10.3yr','p10.5yr')){ # maximum of 3 or 5 year values
      if (is.null(dim(data))){
        out[[att]] = extractor.summaryP10(func=attCalcInfo[[att]]$func,
                                          data=data,
                                          indx=attCalcInfo[[att]]$indx,
                                          attArgs=attCalcInfo[[att]]$attArgs)
      } else {
        out[[att]] = apply(X=data,MARGIN=2,FUN=extractor.summaryP10,
                           func=attCalcInfo[[att]]$func,
                           indx=attCalcInfo[[att]]$indx,
                           attArgs=attCalcInfo[[att]]$attArgs)
      }
    } else if (attCalcInfo[[att]]$opName%in%c('p1.3yr','p1.5yr')){ # maximum of 3 or 5 year values
      if (is.null(dim(data))){
        out[[att]] = extractor.summaryP1(func=attCalcInfo[[att]]$func,
                                          data=data,
                                          indx=attCalcInfo[[att]]$indx,
                                          attArgs=attCalcInfo[[att]]$attArgs)
      } else {
        out[[att]] = apply(X=data,MARGIN=2,FUN=extractor.summaryP1,
                           func=attCalcInfo[[att]]$func,
                           indx=attCalcInfo[[att]]$indx,
                           attArgs=attCalcInfo[[att]]$attArgs)
      }
    }
  }

  return(out)

}

####################################
# Some error handling functions.
# Note: need to fix up use of logfile() - currently just outputting to screen through stop()

invalidSuffixStop = function(funcName,suffix){
  if (is.na(suffix)){
    errMess = paste0("Error: invalid attribute name (must specify suffix for attribute function '",funcName,"')")
  } else {
    errMess = paste0("Error: invalid attribute name (cannot use suffix '",suffix,"' for attribute function '",funcName,"')")
  }
  #  logfile(errMess,file)
  #  logfile("Program terminated",file)
  cat(errMess)
  stop(errMess)
}

invalidStratificationStop = function(strat){
  errMess = paste0("Error: invalid attribute name (stratification '",strat,"' not valid)")
  cat(errMess)
  #  logfile(errMess,file)
  #  logfile("Program terminated",file)
  stop(errMess)
}

invalidOperationStop = function(opName){
  errMess = paste0("Error: invalid attribute name (operation '",opName,"' not valid)")
  cat(errMess)
  #  logfile(errMess,file)
  #  logfile("Program terminated",file)
  stop(errMess)
}

invalidFuncStop = function(func){
  errMess = paste0("Error: invalid attribute name (function '",func,"' does not exist)")
  cat(errMess)
  #  logfile(errMess,file)
  #  logfile("Program terminated",file)
  stop(errMess)
}

####################################
# calculate attribute info based on attribute name

attribute.calculator.setup = function(attSel, # list of evaluated attribute names
                                      datInd #dat indices and properties (e.g. datInd$nyr, datInd$i.yy)
                                      ){

  attCalcInfo = list()

  for (att in attSel){

    # split up attribute name
    chopped=strsplit(x = att,split="_")[[1]]

    # variable name
    varName = chopped[1]
    # stratification index name
    indexName = chopped[2]
    # long function name (including parameters)
    funcNameLong = chopped[3]
    # operator name
    opName = NULL
    if (length(chopped)>3){opName=chopped[4]}

    # calculate selected data indices
    indx = calcStratIndex(indexName,opName,datInd)

    # calculate function names and arguments
    o = calcFuncNamesAndArgs(funcNameLong,datInd)

    attCalcInfo[[att]] = list(func=o$func,attArgs=o$attArgs,indx=indx,opName=opName)

  }

  return(attCalcInfo)

}

####################################
# Calculate function names and arguments

calcFuncNamesAndArgs = function(funcNameLong, # long function name (including parameters)
                                datInd # dat indices and properties (e.g. datInd$nyr, datInd$i.yy)
                                ){

  # functions that require threshold arguments
  funcsWithThresh = c('nWet','dyWet','maxDSD','maxWSD','avgWSD','avgDSD')

  attArgs = NULL

  suffix = NULL
  # select index for funcsWithThresh
  iFuncWithThresh = which(startsWith(funcNameLong,funcsWithThresh))
  if (funcNameLong %in% funcsWithThresh){ # case where no additional params are specified as suffixes
    funcName = funcNameLong
    attArgs=list(threshold=0.) # set default threshold to zero since no additional params specified
  } else if (length(iFuncWithThresh)>0){ # case where suffix is specified
    funcName = funcsWithThresh[iFuncWithThresh]
    suffix = strsplit(funcNameLong,funcName)[[1]][2] # suffix
    if (substring(suffix,1,1)=='T'){ # for case where suffix starts with T we read off threshold
      threshold = as.numeric(substring(suffix,2))
      if (is.na(threshold)){invalidSuffixStop(funcName=funcName,suffix=suffix)} # stop if threshold not numeric
      attArgs = list(threshold=threshold)
    } else {invalidSuffixStop(funcName=funcName,suffix=suffix)} # stop if suffix doesn't start with T

    # seasonality ratios
  } else if (funcNameLong=='seasRatio'){ # note seasRatio not setup to work with monthly/seasonal stratification or with "_m" for mean annual
    funcName = 'seasRatio' # seasonality ratio from foreSIGHT 1.0 (wet season = MAM+JJA)
    phi = (1/12)*2*pi + pi/2.
    attArgs=list(indexWet=c(datInd$i.ss[[3]],datInd$i.ss[[4]]),
                 indexDry=c(datInd$i.ss[[1]],datInd$i.ss[[2]]),
                 phi=phi)
  } else if (startsWith(funcNameLong,'seasRatio')){
    funcName = 'seasRatio'
    suffix = strsplit(funcNameLong,funcName)[[1]][2]
    if (substring(suffix,1,4)=='Mwet'){
      wet1 = match(substring(suffix,5,7),month.abb) # wet season start month
      wet2 = match(substring(suffix,8,10),month.abb) # wet season end month
      if (!is.integer(wet1)|wet1<1|wet1>12|!is.integer(wet2)|wet2<1|wet2>12){invalidSuffixStop(funcName=funcName,suffix=suffix)}
      if (wet1<wet2){
        mWet = seq(wet1,wet2)
      } else { # handle case where wet season ends in next year
        mWet = (seq(wet1,wet2+12)-1)%%12 + 1
      }
      mAll = 1:12
      # calculate dry months
      mDry = mAll[!(mAll%in%mWet)]
      # use middle of dry/wet season to calculate phase of harmonic used in seasonla scaling
      if(max(diff(mDry))==1){
        midDry = stats::median(mDry)-0.5
        monBottom = midDry/12
      } else if (max(diff(mWet))==1) {
        midWet = stats::median(mWet)-0.5
        monBottom = (midWet - 6)/12
      }
      phi = monBottom*2*pi - pi/2
      attArgs=list(indexWet=unlist(datInd$i.mm[mWet]),indexDry=unlist(datInd$i.mm[mDry]),phi=phi)
    } else {
      invalidSuffixStop(funcName=funcName,suffix=suffix)
    }

    # quantile ranges
  } else if (funcNameLong=='rng'){
    funcName = 'rng'
    attArgs=list(lim=0.9) # default limit in 90% (i.e. 5-95%)
  } else if (startsWith(funcNameLong,'rng')){
    funcName = 'rng'
    suffix = strsplit(funcNameLong,funcName)[[1]][2]
    lim = as.numeric(suffix)/100.
    if (!is.na(lim)){
      attArgs=list(lim=lim)
    } else {
      invalidSuffixStop(funcName=funcName,suffix=suffix)
    }

    # percentiles
  } else if (substring(funcNameLong,1,1)=='P'){
    funcName = 'P'
    suffix = strsplit(funcNameLong,funcName)[[1]][2]
    p = as.numeric(suffix)
    if (!is.na(p)){
      attArgs=list(quant=0.01*p) # convert percentile to quantile
    } else {
      invalidSuffixStop(funcName=funcName,suffix=suffix)
    }

    # num days above threshold
  } else if (substring(funcNameLong,1,1)=='R'){
    funcName = 'R'
    suffix = strsplit(funcNameLong,funcName)[[1]][2]
    t = as.numeric(suffix)
    if (!is.na(t)){
      attArgs=list(threshold=t)
    } else {
      invalidSuffixStop(funcName=funcName,suffix=suffix)
    }

  } else if (funcNameLong%in%c('wettest6monPeakDay','wettest6monSeasRatio')){
    funcName = funcNameLong
    doy = datInd$jj
    keepMat = calc_keepMat(doy)
    attArgs = list(doy=doy,keepMat=keepMat)
  } else {
    funcName = funcNameLong
  }

  if (!(funcName%in%attributeFuncs())){invalidFuncStop(func=funcName)}
  func = get(paste0('func_',funcName))

  return(list(func=func,attArgs=attArgs,funcName=funcName,suffix=suffix))

}

####################################
# calculate stratification index

calcStratIndex = function(indexName,opName,datInd){

  # abbreviate for season names
  season.abb = c('SON','DJF','MAM','JJA')
  month.str.abb <- c("JFMAMJJASONDJFMAMJJASOND") #2 year month abbreviation to allow for wrap around months
  month_number <- c(1:12,1:12)   #month.str.abb as month numbers

  stratIndx = NULL
  if (indexName=='ann'){ # this uses all data
    stratIndx = 1:datInd$ndays
  } else if (indexName %in% month.abb){ # this only uses data from given month
    mSel = match(indexName,month.abb)
    stratIndx = datInd$i.mm[[mSel]]
  } else if (indexName %in% season.abb){ # this only uses data from given season
    sSel = match(indexName,season.abb)
    stratIndx = datInd$i.ss[[sSel]]
  } else if (regexpr(indexName,month.str.abb)[1]!=-1 & nchar(indexName) > 1 & nchar(indexName) < 12) { #check if string is sequential months and are greater than 1 month and less than 12
    indexName_posit <- regexpr(indexName,month.str.abb)   #find position of month string in 2 year abbreviation
    mSel <- sort(month_number[seq(indexName_posit[1],length=attr(indexName_posit,"match.length"),by=1)])  #match month letter index to number index and sort months sequentially
    for (n in 1:nchar(indexName)){
      stratIndx = c(stratIndx,datInd$i.mm[[mSel[n]]])
    }
  } else {
    invalidStratificationStop(strat=indexName)
  }

  if (is.null(opName)){
    indx = stratIndx
  } else { # here we calculate stratification for each year (later used to calculate mean/max values over all years)
    yrIndx = list()
    if (opName%in%c('m','sd','cor','dwellTime','range90','corSOI')){
      for (y in 1:length(datInd$i.yy)){
        yrIndx[[y]] = intersect(datInd$i.yy[[y]],stratIndx)
      }
    } else if ((opName=='max3yr')|(opName=='min3yr')|(opName=='p10.3yr')|(opName=='p1.3yr')){
      for (y in 1:length(datInd$i.3yy)){
        yrIndx[[y]] = intersect(datInd$i.3yy[[y]],stratIndx)
      }
    } else if ((opName=='max5yr')|(opName=='min5yr')|(opName=='p10.5yr')|(opName=='p1.5yr')){
      for (y in 1:length(datInd$i.5yy)){
        yrIndx[[y]] = intersect(datInd$i.5yy[[y]],stratIndx)
      }
    } else if (opName=='m10yrBlock'){ # note this is binned average, not moving average (unlike max5yr)
      yrIndx = list()
      for (y in 1:length(datInd$i.10yyBlock)){
        yrIndx[[y]] = intersect(datInd$i.10yyBlock[[y]],stratIndx)
      }
    } else {
      invalidOperationStop(opName=opName)
    }
    indx = yrIndx
  }

  return(indx)

}

####################################
#ATTRIBUTE AUX INFO (determine attribute type and if approved combo with model used)
attribute.info.check<-function(attSel=NULL,  # vector of selected attributes (strings)
                               attPrim=NULL,
                               lambda.mult=NULL
                              #simVar=NULL    # vector of variables simulated using models e.g. c("P","Temp")
                              # modelTag=NULL # model selected
){
  nAtt=length(attSel) # no. of attributes nominated
  attInfo=list()      #create blank list for storage

  #attribute name chopper function
  attInfo$varType=vapply(attSel,FUN = get.attribute.varType,FUN.VALUE=character(1),USE.NAMES = FALSE) #drop use of names as comes ordered anyway

  #ASSIGN TARGET TYPE (IF P USE "FRAC", IF T USE "DIFF")
  attInfo$targetType=vapply(attInfo$varType,FUN=get.target.type,FUN.VALUE=character(1),USE.NAMES=FALSE)

  #FIND WHICH ARE PRIMARY
  if(is.null(attPrim)){
    attInfo$primType=rep(FALSE,nAtt)
    attInfo$primMult=rep(0,nAtt)
  }else{
    get.ind<-function(x,y){which(x == y)}     # quick function to find which are primary attributes
    primInd=vapply(attPrim,FUN=get.ind,FUN.VALUE=numeric(1),x=attSel,USE.NAMES = FALSE)  #Indices of primary attributes
    attInfo$primType=rep(FALSE,nAtt)
    attInfo$primType[primInd]=TRUE   #mark primary ones 'TRUE'

    # Get lambda.mult (corresponding to attPrim) in the correct locations wrt attSel
    # Infomation in lamda.mult will be transferred to primMult & used in the penalty score calculation for the objective function
    # Note: attPrim is a member of attSel - but the order need not be exact, so assignment has to be done inside the loop
    attInfo$primMult=rep(0,nAtt)
    for(i in 1:length(attPrim)){
      # Locate attPrim
      indPrim <- which(attSel == attPrim[i])
      # Assign lambda.mult corresponding to attPrim
      attInfo$primMult[indPrim] <- lambda.mult[i]
    }
  }

  #CHECK FOR INVALID MODEL CHOICE - returns a logical for each attribute
  # attInfo$modelInvalid=vapply(attSel,FUN=check.attribute.model.combo,FUN.VALUE=logical(1),USE.NAMES=FALSE,modelTag=modelTag)

  return(attInfo)
}

get.att.ind<-function(attInfo=NULL,
                      simVar=NULL
){
  #DETERMINE WHICH ATTRIBUTE RELATES TO WHICH SIMULATOR
  attInd=list()
  if(simVar[1] != "All"){                    # ONLY DO IF STOCHASTIC GENERATION IS SELECTED (not simple scaling)
    for(i in 1:length(simVar)){
      attInd[[simVar[i]]]= which(attInfo$varType==simVar[i])
    }
  }
  return(attInd)
}



update.att.Info<-function(attInfo=NULL,
                          attInd=NULL,
                          modelTag=NULL,
                          simVar=NULL
){
  #divide up attInfo to different models
    for(i in 1:length(modelTag)){
      attInfo[[modelTag[i]]]$varType=attInfo$varType[attInd[[simVar[i]]]]
      attInfo[[modelTag[i]]]$targetType=attInfo$targetType[attInd[[simVar[i]]]]
      attInfo[[modelTag[i]]]$primType=attInfo$primType[attInd[[simVar[i]]]]
      attInfo[[modelTag[i]]]$primMult=attInfo$primMult[attInd[[simVar[i]]]]
    }
  return(attInfo)
}

#GETS VARTYPE BY READING FIRST ELEMENT OF ATTRIBUTE STRING
get.attribute.varType<-function(attrib=NULL, # attribute name
                                 sep="_"){
  varType=strsplit(x = attrib,split=sep)[[1]][1]
  return(varType)
}

#get.attribute.varType(attrib=attSel[1], sep="_")

#TARGET TYPE CLASSIFIER
get.target.type<-function(varType=NULL){
  if(varType == "P"){
    targetType="frac"
  }else{
      if(varType == "Temp"){
        targetType="diff"
      }else{
        targetType="frac"
      }
  }
  return(targetType)
}

#######################################################################################

tagBlender<-function(attLab=NULL
){

  # split up attribute name
  chopped=strsplit(x = attLab,split="_")[[1]]

  # variable name
  varName = chopped[1]
  # stratification index name
  indexName = chopped[2]
  # long function name (including parameters)
  funcNameLong = chopped[3]
  # operator name
  opName = NULL
  if (length(chopped)>3){opName=chopped[4]}

  #variable type
  if(varName== "P"){
    vtype="rainfall"
  }else if(varName== "Temp"){
    vtype="temperature"
  }else if(varName== "PET"){
    vtype="PET"
  }else if(varName=="Radn"){
    vtype="Radiation"
  }

  #stratification type
  month.str.abb <- c("JFMAMJJASONDJFMAMJJASOND") #2 year month abbreviation to allow for wrap around months
  month_number <- c(1:12,1:12)   #month.str.abb as month numbers
  if(indexName== "ann"){
    atype="annual"
  }else if(indexName== "DJF"){
    atype="DJF"
  }else if(indexName== "MAM"){
    atype="MAM"
  }else if(indexName== "JJA"){
    atype="JJA"
  }else if(indexName== "SON"){
    atype="SON"
  }else if(indexName== "Jan"){
    atype="Jan"
  }else if(indexName== "Feb"){
    atype="Feb"
  }else if(indexName== "Mar"){
    atype="Mar"
  }else if(indexName== "Apr"){
    atype="Apr"
  }else if(indexName== "May"){
    atype="May"
  }else if(indexName== "Jun"){
    atype="Jun"
  }else if(indexName== "Jul"){
    atype="Jul"
  }else if(indexName== "Aug"){
    atype="Aug"
  }else if(indexName== "Sep"){
    atype="Sep"
  }else if(indexName== "Oct"){
    atype="Oct"
  }else if(indexName== "Nov"){
    atype="Nov"
  }else if(indexName== "Dec"){
    atype="Dec"
  } else if (regexpr(indexName,month.str.abb)[1]!=-1 & nchar(indexName) > 1 & nchar(indexName) < 12) {
    atype=indexName
  } else {
    cat(paste0('invalid attribute: cannot use ',indexName,' stratification'))
    return(invisible())
  }

  # use calcFuncNamesAndArgs() to calculate parameter values from long function name
  o = calcFuncNamesAndArgs(funcNameLong = funcNameLong,datInd = NULL)
  if(funcNameLong== "tot"){
    mtype="total"
  } else if(funcNameLong== "avg"){
    mtype="average"
  } else if (startsWith(funcNameLong,'seasRatio')){
    if (is.null(o$suffix)){
      mtype='ratio of wet to dry season totals'
    } else {
      wetStart = substring(o$suffix,5,7)
      wetEnd = substring(o$suffix,8,10)
      mtype=paste0('ratio of wet (',wetStart,'-',wetEnd,') to dry season totals')
    }
    if(indexName=='ann'){
      atype = NULL
    } else {
      errMess = paste0('invalid attribute: cannot compute seasRatio for ',indexName,' stratification\n')
      cat(errMess)
      return(invisible())
    }
  } else if (substring(funcNameLong,1,1)=='P'){
    if (is.null(o$suffix)){
      errMess = 'invalid attribute: P attribute requires specification of percentile\n'
      cat(errMess)
      return(invisible())
    } else {
      p=o$suffix
      mtype=paste0(p,'th percentile')
    }
    if(indexName== "ann"){atype=NULL}
  } else if (startsWith(funcNameLong,'nWet')){
    if (is.null(o$suffix)){
      mtype="no. wet days"
    } else {
      thresh = o$suffix
      mtype=paste0('no. wet days (above ',thresh,')')
    }
  } else if (startsWith(funcNameLong,'maxDSD')){
    if (is.null(o$suffix)){
      mtype="max dryspell duration"
    } else {
      thresh = o$suffix
      mtype=paste0('max dryspell duration (below ',thresh,')')
    }
  } else if (startsWith(funcNameLong,'maxWSD')){
    if (is.null(o$suffix)){
      mtype="max wetspell duration"
    } else {
      thresh = o$suffix
      mtype=paste0('max wetspell duration (above ',thresh,')')
    }
  } else if (startsWith(funcNameLong,'avgDSD')){
    if (is.null(o$suffix)){
      mtype="average dryspell duration"
    } else {
      thresh = o$suffix
      mtype=paste0('average dryspell duration (below ',thresh,')')
    }
  } else if (startsWith(funcNameLong,'avgWSD')){
    if (is.null(o$suffix)){
      mtype="average wetspell duration"
    } else {
      thresh = o$suffix
      mtype=paste0('average wetspell duration (above ',thresh,')')
    }
  } else if (startsWith(funcNameLong,'dyWet')){
    if (is.null(o$suffix)){
      mtype="wet day amount"
    } else {
      thresh = o$suffix
      mtype=paste0('wet day amount (above ',thresh,')')
    }
  }else if(funcNameLong== "GSL"){
    mtype="growing season length"
  }else if(funcNameLong== "CSL"){
    mtype="cold season length"
  }else if(funcNameLong== "F0"){
    if(varName!='T'){
      errMess = 'invalid attribute: can only compute frost days for T\n'
      cat(errMess)
      return(invisible())
    }
    mtype="no. frost days"
  } else if (substring(funcNameLong,1,1)=='R'){
    if (is.null(o$suffix)){
      errMess = 'invalid attribute: R attribute requires specification of threshold\n'
      cat(errMess)
      return(invisible())
    } else {
      t=o$suffix
      mtype=paste0('no. days above ',t)
    }
  } else if (startsWith(funcNameLong,'rng')){
    if (is.null(o$attArgs$lim)){
      errMess = 'invalid attribute: rng attribute requires attrtibe argument for lim\n'
      cat(errMess)
      return(invisible())
    } else {
      lim = 100*as.numeric(o$attArgs$lim)
      mtype=paste0(lim,'% range')
    }
  } else {
    errMess = paste0('invalid attribute: built-in function not available for ',funcNameLong)
    cat(errMess)
    return(invisible())
  }

  #statType
  if(is.null(opName)){
    stype = ''
  } else if (opName=='m'){
    stype="Mean"
  } else if (opName=='sd'){
    stype="Sdev"
  } else if (opName=='min5yr'){
    stype="Min 5yr total"
  } else if (opName=='dwellTime'){
    stype="dwell time"
  } else if (opName=='range90'){
    stype="90% range"
  } else {
    stype=opName
  }

  #stitch togther
  phrase=paste(stype,atype,mtype,vtype)
  return(phrase)
}


