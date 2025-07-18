rm(list=ls())

foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT/'
devtools::load_all(foreSIGHTDir)

dirname = 'C:/Users/a1065639/Box/2025_DEW_foreSIGHT/Data/SA_climate_ready/'

#site = '23343'
site = '23090'
region = 'amlr27'
#modelList = c('access10','bcc.csm11m','gfdl.esm2m','ipsl.cm5blr','miroc5','mri.cgcm3')
modelList = c('canesm2','cnrm.cm5','gfdl.esm2m','ipsl.cm5blr','miroc5')#,'mri.cgcm3')
#period = 'his'
#scenList = c('his','r45','r85')
scenList = c('his','r85')

attSel = c('P_day_all_tot_m','P_day_all_seasRatioMarAug',
           'P_day_all_P99','P_day_all_avgDSD','P_day_all_nWet_m',
           'P_day_all_tot_cv','P_day_all_tot_dwellTime',
           'Temp_day_all_avg')
nReps = 3

year.start.his = 1976
year.end.his = 2000

year.start.fut = 2076
year.end.fut = 2100

systemArgs = systemArgs <- list(roofArea = 205, 
                                nPeople = 1, 
                                tankVol = 2400, 
                                firstFlush = 2.0, 
                                write.file = FALSE)

metrics <- c("average daily deficit (L)", "reliability (fraction)", "volumetric reliability (fraction)")

df.atts <- data.frame(matrix(ncol = (3+length(attSel)+length(metrics)), nrow = 0))

# df.atts.bias <- data.frame(matrix(ncol = (2+length(attSel)), nrow = 0))

# fname = 'C:/Users/a1065639/Work/foreSIGHT/testing_2.0/barossa_clim_ref.RData'
# load(fname)
# clim_ref$P = apply(clim_ref$P,1,mean)
# attsClim = calculateAttributes(clim_ref,attSel)

for (model in modelList){
  
  for (r in 1:nReps){

    print(r)
    rep = sprintf("%03d", r)
    
    for (s in 1:length(scenList)){
      
      scen = scenList[s]
      
      fname = paste0(dirname,site,'/',model,'/',region,'.',model,'.',scen,'.',site,'.',rep,'.txt')
      
      dat = read.table(fname)
      
      df.dates <- data.frame(
        day = as.integer(dat[,3]),
        month = as.integer(dat[,2]),
        year = as.integer(dat[,1])
      )
      
      dates <- as.Date(with(df.dates, paste(year, month, day, sep = "-")))
      
      rain = as.numeric(dat[,5])
      temp.max = as.numeric(dat[,6])
      temp.min = as.numeric(dat[,7])
      temp = 0.5*(temp.min+temp.max)
      
      if (scen=='his'){
        keep = which((df.dates$year>=year.start.his)&(df.dates$year<=year.end.his))
      } else {
        keep = which((df.dates$year>=year.start.fut)&(df.dates$year<=year.end.fut))
      }
      
      clim = list(times=as.POSIXct(dates[keep],tz='UTC'),
                  P=rain[keep],Temp=temp[keep])
      
      metricVals = unlist(tankWrapper(data=clim,systemArgs=systemArgs,metrics=metrics))
      
      attVals = calculateAttributes(clim,attSel)
      
      if (scen=='his'){
        attValsHis = attVals
        # attValsBias = attVals/attsClim-1
        # new.row = c(model,rep,attValsBias)
        # df.atts.bias = rbind(df.atts.bias,new.row)
      } else {
        attValsFac = attVals/attValsHis
        new.row = c(attValsFac,scen,model,rep,metricVals)
        df.atts = rbind(df.atts,new.row)
      }
      
    }
  }
  
}

colnames(df.atts) <- c(attSel,'scen','model','rep',metrics)

numericCol = c( 1:length(attSel) , (length(attSel)+4):ncol(df.atts))

df.atts[numericCol] = as.numeric(as.matrix(df.atts[,numericCol]))

df.atts = subset(df.atts,select=-rep)
df.atts = subset(df.atts,select=-scen)
df.colnames = colnames(df.atts) 
colnames(df.atts)[df.colnames=='model'] = 'Name'

egClimDataDM = df.atts

save(file='data/egClimDataDM.rda',egClimDataDM)

# climData = data.frame(P_day_all_tot_m=as.numeric(df.atts$P_day_all_tot_m),
#                       P_day_all_seasRatioMarAug=as.numeric(df.atts$P_day_all_seasRatioMarAug),
#                       Temp_day_all_avg=as.numeric(df.atts$Tmax_day_all_avg),
#                       Name=df.atts$model)

# P_ann_tot_m P_ann_seasRatio P_ann_nWet_m Temp_ann_avg_m    Name Avg. Deficit
# 1   0.8978995       1.1138472    0.5118911      0.5810644 Model-C           26
# 2   0.9702083       0.9211320    0.5561901      0.5523804 Model-C           31
