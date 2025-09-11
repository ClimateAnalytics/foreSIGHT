########################

foreSIGHTDir = 'C:/Users/a1065639/Work/foreSIGHT/'
droughRiskDir = 'C:/Users/a1065639/Work/DroughtRisk/'
dataDir.Loxton = 'C:/Users/a1065639/Box/2025_DEW_foreSIGHT/Data/'
dataDir.Barossa = 'C:/Users/a1065639/Box/2021 CRAFT Barossa project/External datasets/Hydroclimate (instrumental)/'

# foreSIGHTDir = '../'
# droughRiskDir = '/hpcfs/users/a1065639/git/DroughtRisk_Feb_2024_paper_shared/DroughtRisk/'
# dataDir.Loxton = '../../Data/'
# dataDir.Barossa = '../../Data/'

devtools::load_all(foreSIGHTDir)
devtools::load_all(droughRiskDir)

runDirname = paste0(foreSIGHTDir,'testing_2.0/')

source('testing_2.0/load_camels.R')