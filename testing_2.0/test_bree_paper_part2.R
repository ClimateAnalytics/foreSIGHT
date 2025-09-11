rm(list=ls())

library(foreSIGHT)

setwd("C:/Users/a1065639/Work/foreSIGHT/data/scripts/")

print('loading data')
load('bree_test.RData')

print('getting sim summary')
simSummary = getSimSummary(simTest)

systemArgsA <- list(roofArea = 205,                         # roof area in m2
                    nPeople = 1,                            # number of people using water
                    tankVol = 2400,                         # tank volume in L
                    firstFlush = 2.0,                       # first flush volume
                    write.file = F)
metrics = c("average daily deficit (L)", "reliability (fraction)") # performance metrics chosen

print('calcuating system performance')
systemA.perf <- runSystemModel(sim=simTest,
                               systemModel=tankWrapper,
                               systemArgs=systemArgsA,
                               metrics = metrics)

plotPerformanceSpace(performance = systemA.perf[2], sim = simSummary)


#data("egSimSummary")
#data("egSimPerformance") 
data("egClimData")
# plot performance space using the metric average deficit
threshAvgDeficit <- 28
labelAvgDeficit <- "Max Avg. Deficit (28L)"
defColLim <- c(17, 33)
topReps <- 7
#plotPerformanceSpace(egSimPerformance[1], egSimSummary, attX = "P_ann_seasRatio_m", attY = "P_ann_tot_m",  topReps = topReps, colLim = defColLim, perfThresh = threshAvgDeficit, perfThreshLabel = labelAvgDeficit, climData = egClimData)
print('plotting system performance')
plotPerformanceSpace(systemA.perf[1], simSummary, attX = "P_day_all_seasRatio", attY = "P_day_all_tot_m",  topReps = topReps, colLim = defColLim, perfThresh = threshAvgDeficit, perfThreshLabel = labelAvgDeficit, climData = egClimData)


pause

# plot performance space using the metric reliability
threshRel <- 0.81
labelRel <- "Min Reliability (0.81)"
relColLim <- c(0.76, 0.88)
topReps <- 7
plotPerformanceSpace(egSimPerformance[2], egSimSummary, attX = "P_ann_seasRatio_m", attY = "P_ann_tot_m",  topReps = topReps, colLim = relColLim, perfThresh = threshRel, perfThreshLabel = labelRel, climData = egClimData)
# This is a maximum thershold
threshAvgDeficit <- 28
# This is a minimum threshold
threshRel <- 0.81
perfThreshMax <- c(threshAvgDeficit, NA)
perfThreshMin <- c(NA, threshRel)
topReps <- 7
plotPerformanceSpaceMulti(egSimPerformance, egSimSummary, attX = "P_ann_seasRatio_m", attY = "P_ann_tot_m", topReps = topReps, perfThreshMin = perfThreshMin, perfThreshMax = perfThreshMax, climData = egClimData, col = viridisLite::inferno(6, direction = -1))

