rm(list=ls())

devtools::load_all()


########################################

modelSelection = list()
modelSelection$modelType = list()
modelSelection$modelType$P = "distScaling"
modelSelection$modelParameterVariation = list()
modelSelection$modelParameterVariation$P = "annual"
modelSelectionJSON = jsonlite::toJSON(modelSelection, pretty = TRUE, auto_unbox = TRUE)
controlFile = paste0(tempdir(), "\\eg_controlFile.json")
write(modelSelectionJSON, file = controlFile)

########################################

attPerturb = c("P_day_all_P99")
attHold = c("P_day_all_tot_m")

attPerturbType = "regGrid"
attPerturbSamp = c(1)
attPerturbMin = c(1.5)
attPerturbMax = c(1.5)

# create the exposure space
expSpace = createExpSpace(attPerturb = attPerturb,
                          attPerturbSamp = attPerturbSamp,
                          attPerturbMin = attPerturbMin,
                          attPerturbMax = attPerturbMax,
                          attPerturbType = attPerturbType,
                          attHold = attHold)

########################################

sim = generateScenarios(reference = tank_obs,
                              expSpace = expSpace,
                              controlFile = controlFile,
                              seedID = 1)

plotScenarios(sim)
