rm(list = ls())

source('replaceStrFile.R')

catchmentList = c('Loxton','Barossa.1','410057','A5030502','419005','410730')

attPerturbList = c("P_day_all_tot", "P_day_all_xP99overPave", "P_day_all_avgDSD", "PET_day_all_avg")

startYr = 1976
endYr = 2005

numReplicates = 10
cores = 25

submitFile = 'template_submit_OAT_25'

for (catchment in catchmentList){
  print(catchment)
  for (attPerturb in attPerturbList){
    print(attPerturb)

    newSubmitFile = paste0(submitFile,'_',catchment,'_',attPerturb,'.tmp',sep='')
    file.copy(from=submitFile,to=newSubmitFile)
    replaceStrFile(newSubmitFile,'catchment',catchment)
    replaceStrFile(newSubmitFile,'attPerturb',attPerturb)

    command = paste('sbatch ', newSubmitFile, sep='')
    print(command)
    system(command)

    file.remove(newSubmitFile)

  }
}






