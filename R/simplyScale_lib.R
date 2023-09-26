#######################################
##   SIMPLE SCALING FUNCTIONS     ##
#######################################

# CONTAINS
  #function for scaling over a period

simple.scaling<-function(target=NULL,          #extracted from matrix output of exposure space maker. annual format is a vector, passed from wrapper.
                         targetType=NULL,      #diff or frac
                         data=NULL,            #data series stripped of dates
                         varType=NULL,         #names of data series
                         period=NULL,          #number of separations for scaling e.g. annual all at once, seasonal in four stages.
                         i.pp=NULL            #index to control what time step entries are changed at a time
) {

   temp = list()
   nLoc = list()
   for (var in varType){
     dimVar = dim(data[[var]])
     nLoc[[var]] = dimVar[2]
     temp[[var]] = matrix(NA,nrow=dimVar[1],ncol=dimVar[2])
     colnames(temp[[var]])=colnames(data[[var]])
    }

   for (p in 1:period) {           #for annual this is equal to 1
    #select target subset if doing seasonal
    #******need method of dealing with different targets
     for (j in 1:length(varType)) {              #loop for each variable, to check its scale type.
       var = varType[j]
       switch(targetType[j],
             "frac" = {temp[[var]][i.pp[[p]],1:nLoc[[var]]]=as.matrix(data[[var]][i.pp[[p]],1:nLoc[[var]]])*target[j]},
             "diff" = {temp[[var]][i.pp[[p]],1:nLoc[[var]]]=as.matrix(data[[var]][i.pp[[p]],1:nLoc[[var]]])+target[j]},
             -99.00)
    }
  }
  return(temp)
}

###################################################################

simple.scaling.seasonal<-function(target_total_fac=NULL, # multiplication factor for tot/ave
                                  target_seas_fac=NULL,  # multiplication factor for seas ratio
                                  varType=NULL,          # variable name (P, PET)
                                  data=NULL,             # data time series
                                  i.T=NULL,              # indices for all data
                                  i.S1=NULL,             # indices for dry period
                                  i.S2=NULL,             # indices for wet period
                                  phi=NULL               # phase of harmonic
) {

  if (length(varType)>1){
    stop("Incorrect use of argument varType in simple.scaling.seasonal()")
  } else if (!(varType%in%c('P','PET'))) {
    stop("Seasonal scaling only setup for P and PET")
  } else {
    var = varType
  }

  temp = list()
  nLoc = list()
  dimVar = dim(data[[var]])
  nLoc[[var]] = dimVar[2]
  temp[[var]] = matrix(NA,nrow=dimVar[1],ncol=dimVar[2])
  colnames(temp[[var]])=colnames(data[[var]])

  date = as.Date(paste0(data$year,'/',data$month,'/',data$day))
  dateOneYear = date[1:365]
  doy = as.integer(format(date,'%j'))

#  alpha = target_total_fac
#  beta = target_seas_fac

#  phiLoc = phi

  for (l in 1:nLoc[[var]]){

    # DM: This code allows for phase to be computed from observed data - currenty not enabled
    # if (is.null(phi)){
    #   clim = calc_ClimDaily_dayOfYearWindow(obs=data[[var]][i.T,l],dateObs = date,dateClim=dateOneYear,inc = 14)
    #   meanClim = apply(clim,1,mean)
    #   o = HarmonicRegression::harmonic.regression(inputts = meanClim, inputtime = seq(1,365),Tau = 365)
    #   phiLoc = o$pars$phi
    # } else {
    #  phiLoc = phi
    # }

    P = data[[var]][i.T,l]
    H = sin(phi+2*pi*doy/365)
    PH = P*H

    I1_T = sum(P[i.T])
    I2_T = sum(PH[i.T])
    I1_S1 = sum(P[i.S1])
    I2_S1 = sum(PH[i.S1])
    I1_S2 = sum(P[i.S2])
    I2_S2 = sum(PH[i.S2])

    RHS = vector(length=2)
    coeff = matrix(nrow=2,ncol=2)

    coeff[1,1] = I1_T
    coeff[1,2] = I2_T
    RHS[1] = target_total_fac*I1_T

    coeff[2,1] = (target_seas_fac-1)*I1_S2*I1_S1
    coeff[2,2] = target_seas_fac*I1_S2*I2_S1 - I1_S1*I2_S2
    RHS[2] = 0

    o = solve(a=coeff,b=RHS)

    A = o[1]
    B = o[2]

    P_star = P*(A+B*H)

#    alpha_est = mean(P_star)/mean(P)

#    SR = mean(P[i.S2])/mean(P[i.S1])
#    SR_star = mean(P_star[i.S2])/mean(P_star[i.S1])

#    beta_est = SR_star / SR

    temp[[var]][i.T,l] = P_star

    if(any(P_star<0.)){stop()}

  }
  return(temp)
}

#############################
# climatology based on window about day of year (see McInerney et al, WRR, 2020)
calc_ClimDaily_dayOfYearWindow = function(obs,dateObs,dateClim,inc){

  dObs = as.integer(format(dateObs,'%j'))
  leapYear = as.integer(is.leapyear(as.integer(format(dateObs,'%Y'))))
  nYear = round(length(dateObs)/365.25)

  day_clim = matrix(nrow=366,ncol=(2*inc+1)*nYear)
  for (d in 1:365){
    in_window = (abs(dObs-d)<=inc) |
      (dObs-d)>=(365+leapYear-inc) |
      (dObs-d+365+leapYear)<=inc
    day_clim[d,] = obs[in_window]
  }
  day_clim[366,] = day_clim[365,]

  nDaysClim = length(dateClim)
  clim = matrix(nrow=nDaysClim,ncol=(2*inc+1)*nYear)
  for (n in 1:nDaysClim){
    d = as.integer(format(dateClim[n],'%j') )
    clim[n,] = day_clim[d,]
  }

  return(clim)
}

# #############################

is.leapyear=function(year){
  #http://en.wikipedia.org/wiki/Leap_year
  return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))
}

###################################################################

#TEST
#simple.scaling(target=target,targetType=targetType, data=mat,varType=varType,period=1,index=i.pp )

#----------------------------------------------
# try to imitate format seen in generic generation functions in WGEN_lib.R
