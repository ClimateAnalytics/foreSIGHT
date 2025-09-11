rm(list=ls())

library(dplyr)
library(lubridate)

#########################################################

calculateBreakSeason = function(dates,flow,threshold,minMon=4,maxMon=12){
  
  df <- data.frame(date = dates, value = flow)
  
  # Function to get first day above threshold and days after start
  df_result <- df %>%
    mutate(
      year = year(date),
      month = month(date)
    ) %>%
    # Keep only months between min_month and max_month (inclusive)
    filter(month >= min_month & month <= max_month) %>%
    arrange(date) %>%
    group_by(year) %>%
    mutate(
      cum_sum = cumsum(value),
      period_start = as.Date(paste0(year, "-", sprintf("%02d", min_month), "-01")),
      days_after_start = as.numeric(date - period_start)
    ) %>%
    filter(cum_sum > threshold) %>%
    slice_head(n = 1) %>%
    select(year, first_exceed_date = date, cum_sum, days_after_start) %>%
    ungroup()

  return(df_result$days_after_start)
  
}

#########################################################

calculateNumFlowsYear = function(dates,flow,threshold){
  
  df <- data.frame(date = dates, value = flow)
  
  df_days_above <- df %>%
    mutate(
      year = year(dates),
    ) %>%
    group_by(year) %>%
    summarise(
      days_above_threshold = sum(value > threshold),
      .groups = "drop"
    )
  
  return(df_days_above$days_above_threshold)
  
}

#########################################################

library(zoo)  # for rolling functions

calculateMultiYearTot = function(dates,flow,nYears){
  
  df <- data.frame(date = dates, value = flow)
  
  df_yearly <- df %>%
    mutate(year = year(date)) %>%
    group_by(year) %>%
    summarise(year_total = sum(value), .groups = "drop")
  
  df_yearly <- df_yearly %>%
    mutate(rolling_multiYr_total = rollapply(year_total, width = nYears, FUN = sum, align = "right", fill = NA))
  
  min_total <- min(df_yearly$rolling_multiYr_total, na.rm = TRUE)
  
  df_min <- df_yearly %>%
    filter(rolling_multiYr_total == min_total)
  
  return(df_min$rolling_multiYr_total)
  
}

#########################################################


dirname = 'C:/Users/a1065639/Work/foreSIGHT/testing_2.0/OneDrive_1_18-07-2025/tot/tot/'

flowThreshold = 0.04

# Parameters
min_month <- 4   # April
max_month <- 12  # December

siteNum = 17 # upper jacobs, A5050533
cumFlowThreshold = 95

# siteNum = 18 # upper tanunda
# cumFlowThreshold = 50

nReps = 2
nTars = 2

avgBreakSeason = avgNumFlowsYear = meanAnnualFlow = P99 = P25 = minTot3Year = matrix(nrow=nTars,ncol=nReps)

for (rep in 1:nReps){
  for (tar in 1:nTars){
   
    outFilename = paste0(dirname,'out_Rep',rep,'_Target',tar,'.res.csv')
    sourceOutData = read.csv(outFilename,skip=30,header=F)
    outColNames <-  colnames(read.csv(outFilename,skip = 28,header = T))                  #put data column names into outColNames #skip needs to be based on no col. 
    names(sourceOutData) <- outColNames
    
    dates = as.Date(sourceOutData$Date)
    flow = sourceOutData[,(siteNum+1)]
    
    breakSeason = calculateBreakSeason(dates,flow,cumFlowThreshold,minMon=4,maxMon=12)
    avgBreakSeason[tar,rep] = mean(breakSeason)
    
    NumFlowsYear = calculateNumFlowsYear(dates,flow,flowThreshold)
    avgNumFlowsYear[tar,rep] = mean(NumFlowsYear)
    
    meanAnnualFlow[tar,rep] = mean(flow,na.rm=T)*365.25
    
    P99[tar,rep] = quantile(flow,probs = 0.99,na.rm=T)

    P25[tar,rep] = quantile(flow,probs = 0.25,na.rm=T)
    
    Tot3Year = calculateMultiYearTot(dates,flow,nYears=3)
    minTot3Year[tar,rep] = min(Tot3Year)
    
  }
}

apply(avgBreakSeason,1,median)
apply(avgNumFlowsYear,1,median)
apply(meanAnnualFlow,1,median)
apply(P99,1,median)
apply(P25,1,median)
apply(minTot3Year,1,median)




