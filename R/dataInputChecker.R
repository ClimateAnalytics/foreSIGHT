#########################################
##   CHECKS FOR DATA INPUT TIMESERIES  ##
#########################################

# CONTAINS
  #input checking for data frame

input_check <- function(obs,
                        file,
                        simLengthNyrs=NULL
                        ){            #This input has its dates in first three columns

  # convert obs from df to list. this enables backward-compatability with df inputs.  
  if(is.data.frame(obs)){obs=as.list(obs)}
  # convert data vectors to matrices so that all are matrices. 
  obsVars <- names(obs)[-which(names(obs) %in% c("year", "month", "day"))]
  for (var in obsVars){
    if(is.vector(obs[[var]])){obs[[var]]=as.matrix(obs[[var]],ncol=1)}
  }
  
  nT = length(obs$year)
  keep = seq(1:nT)
    
  #TRUNCATE TO START AND END AT WHOLE YEAR
  first=obs$year[1]
  if(obs$month[1]!=1 && obs$day[1]!=1){
    keep<-keep[obs$year!=first]
  }
  last=obs$year[nT]
  if(obs$month[nT]!=12 && obs$day[nT]!=31){
    keep<-keep[obs$year!=last]
  }
  nT = length(keep)

  for (var in names(obs)){
    if (is.null(dim(obs[[var]]))){
      obs[[var]] = obs[[var]][keep]
    } else {
      obs[[var]][] = obs[[var]][keep,]
    }
    obs[[var]][obs[[var]] <= -99]<-NA
  }
  
  #HOW MANY VARIABLES
  n=length(names(obs))
  
  #ERROR VALUES/MISSING VALUES
  M <- sapply(obs, function(x) sum(is.na(x))) # this might need to be fixed up when we have multiple cats and some missing data
  for (Var in names(M)) {
    if (M[Var]>0) {
      entry <- which(is.na(obs[,Var]))
      dates <- cbind(obs[entry,"year"],obs[entry,"month"],obs[entry,"day"])
      warn(p("Missing entries in ",Var),file)
      cat(apply(dates,1,paste,collapse="-"),sep="\n")
      cat("\n")
      logfile(apply(dates,1,paste,collapse="-"),file)
      
    }
  }
  
  #EXTRACT FIRST DAY AND LAST DAY
  dateS=paste(obs$year[1],str(obs$month[1],2,"0"),str(obs$day[1],2,"0"),sep="-")
  dateF=paste(obs$year[nT],str(obs$month[nT],2,"0"),str(obs$day[nT],2,"0"),sep="-")
  temp=as.numeric(length(date_gen <- seq(as.Date(dateS),as.Date(dateF),by="day")))
  
  
  #CHECK IF SIMLENGTHNYRS> SUPPLIED DATA
  if(!is.null(simLengthNyrs)){
    last=obs$year[1]+simLengthNyrs-1
    obsLast=obs$year[nT]
    if(last<obsLast){ 
      logfile("Error: Simulation length requested is shorter than observed data",file)
      logfile("Program terminated",file)
      stop("Simulation length requested is shorter than observed data")
    }
  }
  
  
  if(sum(M)>0){
    logfile("Error: There are missing data entries in the variables",file)
    logfile("Program terminated",file)
    stop("There are missing data entries in the variables")
    
  } 
  if (temp<3650){
    warn("You have provided less than 10 years of data",file)
  } 
  if (nT!=temp) {
    logfile("Error: There are missing dates from the provided data. Ensure leap days are included",file)
    logfile("Program terminated",file)
    stop("There are missing dates from the provided data. Ensure leap days are included")
  }
  
  return(list(data=obs))
} #END
