rm(list=ls())

# foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT/'
# devtools::load_all(foreSIGHTDir)

library(foreSIGHT)

# directory where SA climate ready data is stored
dirname = 'C:/Users/a1065639/Box/2025_DEW_foreSIGHT/Data/SA_climate_ready/'

# site and region
#site = '23343' # rosedale
site = '23090' # kent town
region = 'amlr27'

# list of GCMs
modelList = c('canesm2','cnrm.cm5','gfdl.esm2m','ipsl.cm5blr','miroc5')

# climate scenarios
scenList = c('his','r85')

# climate attributes to be calculated and saved in dataframe
attSel = c('P_day_all_tot_m',
           'P_day_all_seasRatioMarAug',
           'P_day_all_P99',
           'Temp_day_all_avg')

# number of replicates
nReps = 3

# start and end dates for historical period
year.start.his = 1976
year.end.his = 2000

# start and end dates for future period
year.start.fut = 2076
year.end.fut = 2100

# system model configurations
system = 'A'
#system = 'B'

if (system=='A'){
  systemArgs = systemArgs <- list(roofArea = 205, 
                                  nPeople = 1, 
                                  tankVol = 2400, 
                                  firstFlush = 2.0, 
                                  write.file = FALSE)  
} else if (system=='B'){
  systemArgs <- list(roofArea = 205,
                      nPeople = 1, 
                      tankVol = 2600,
                      firstFlush = 2.0,
                      write.file = FALSE)  
}

# metrics ot be evlauated and saved in dataframe
metrics = c("average daily deficit (L)")
metricNames = c('Avg. Deficit')

# setup dataframe
df.atts <- data.frame(matrix(ncol = (3+length(attSel)+length(metrics)), nrow = 0))

for (model in modelList){
  
  for (r in 1:nReps){

    print(r)
    rep = sprintf("%03d", r)
    
    for (s in 1:length(scenList)){
      
      scen = scenList[s]

      # load ascii data files      
      fname = paste0(dirname,site,'/',model,'/',region,'.',model,'.',scen,'.',site,'.',rep,'.txt')
      dat = read.table(fname)
      
      # read dates
      df.dates <- data.frame(
        day = as.integer(dat[,3]),
        month = as.integer(dat[,2]),
        year = as.integer(dat[,1])
      )
      dates <- as.Date(with(df.dates, paste(year, month, day, sep = "-")))
      
      # read climate data
      rain = as.numeric(dat[,5])
      temp.max = as.numeric(dat[,6])
      temp.min = as.numeric(dat[,7])
      temp = 0.5*(temp.min+temp.max)
      
      # subset data 
      if (scen=='his'){
        keep = which((df.dates$year>=year.start.his)&(df.dates$year<=year.end.his))
      } else {
        keep = which((df.dates$year>=year.start.fut)&(df.dates$year<=year.end.fut))
      }
      
      # setup reference climate data
      clim = list(times=as.POSIXct(dates[keep],tz='UTC'),
                  P=rain[keep],Temp=temp[keep])
      
      # calculate system metrics (and convert to vector)
      metricVals = unlist(tankWrapper(data=clim,systemArgs=systemArgs,metrics=metrics))
      
      # calculate attributes
      attVals = calculateAttributes(clim,attSel)
      
      # compute baseline from historical
      if (scen=='his'){
        attValsHis = attVals
      } else {
        # compute changes in attributes
        attValsFac = attVals/attValsHis
        new.row = c(attValsFac,scen,model,rep,metricVals)
        df.atts = rbind(df.atts,new.row)
      }
      
    }
  }
  
}

# format dataframe
colnames(df.atts) <- c(attSel,'scen','model','rep',metricNames)
numericCol = c( 1:length(attSel) , (length(attSel)+4):ncol(df.atts))
df.atts[numericCol] = as.numeric(as.matrix(df.atts[,numericCol]))
df.atts = subset(df.atts,select=-rep)
df.atts = subset(df.atts,select=-scen)
df.colnames = colnames(df.atts) 
colnames(df.atts)[df.colnames=='model'] = 'Name'

# save dataframe to rda file
if (system=='A'){
  egClimData = df.atts
  save(file='data/egClimData.rda',egClimData)
} else if (system=='B'){
  egClimDataB = df.atts
  save(file='data/egClimDataB.rda',egClimDataB)
}
