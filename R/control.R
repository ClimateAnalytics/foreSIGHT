# spaceInfo=expSpaceSampManager(exSpArgs=exSpArgs,attInfo=attInfo,attSel=attSel,file=file,IOmode=IOmode,nSeed=nSeed,seedID=seedID,arrayID=arrayID)  #allows targetMat to be provided as CSV
# targetMat=spaceInfo$targetMat
# attRot=spaceInfo$attRot      #attRot - used if needed in "OAT" operations
# seedCatalogue=spaceInfo$seedCatalogue   #create list of seeds
# nTarget=dim(targetMat)[1]
#
# #seed generation can go here, making seedCatalogue a matrix, ntarget rows by ndays in dev mode.
# progress("Exposure space sampled OK",file)
# progress(p(nTarget," targets sampled across ", length(attInfo$varType)," attributes"),file)

#-------------------------------------------------------------------------
# Wrapper function that loops over generateScenario
# Kept this function very simple to a loop over targets & replicates - so that it can be opened up for paralellization


#' Produces time series of hydroclimatic variables for an exposure space.
#'
#' \code{generateScenarios} produces time series of hydroclimatic variables using requested climate attributes that correspond to a target exposure space using a reference daily time series as an input.
#' @param reference data.frame or list; contains reference daily climate data. \cr
#'            For single site data, \code{reference} is a data.frame with columns named \emph{year}, \emph{month}, \emph{day}, \emph{*variable_name1*}, \emph{*variable_name2*}.
#'            Note that the first three columns of the data.frame contain the year, month, and day of the data.
#'            The columns have to be named as specified.
#'            For multi-site data, \code{reference} is a list, with elements named \emph{year}, \emph{month}, \emph{day}, \emph{*variable_name1*}, \emph{*variable_name2*}. List format is suitable for both single and multi-site data.
#'            Climate variables are specified as matrices, with columns for each site. \cr
#'            Use \code{viewModels()} to view the valid variable names.
#'            Please refer to data provided with the package that may be loaded using \code{data(tankDat)} and \code{data(barossaDat)} for examples of the expected format of single site and multi-site \code{reference}.
#' @param expSpace a list; created using the function \code{createExpSpace}
#' @param simLengthNyrs a number; a scalar that specifies the length in years of each generated scenario. This argument is used only with stochastic generation.
#' If \code{NULL} (the default), the generated simulation will be as long as \code{reference}.
#' @param numReplicates a number; a scalar that specific the number of stochastic replicates to be generated. The default is 1.
#' @param seedID a number; a scalar that specifies the seed to be used for the first replicate. Subsequent replicates will use seeds incremented by one.
#'               If \code{seedID} is \code{NULL} (which is the default), the function will use a random seed for stochastic time series generation.
#'               The seed used will be specified in the output. This argument is intended for use in cases that aim to reproduce an existing simulation.
#' @param controlFile a string; to specify the model/optimisation options used for simulating time series data. The valid values are:
#' \itemize{
#' \item {\code{NULL}} {: the simulation uses the foreSIGHT default stochastic model settings.}
#' \item {\code{"scaling"}} {: the simulation uses scaling (simple/seasonal) instead of a stochastic model.
#'                           If all attributes in \emph{expSpace} are annual totals/averages, then simple scaling is used.
#'                           If seasonality ratio attributes are also included in \emph{expSpace}, then seasonal scaling is used.}
#' \item {\code{path to a JSON file}} {: the JSON file contains advanced options specify the stochastic model and optimisation inputs.
#'                   These options can be used to change stochastic model types, overwrite default model parameter bounds, change default optimisation arguments, and set penalty attributes to be used in optimisation.
#'                   Please refer to the function \code{writeControlFile} in order to create an \code{controlFile} JSON file.
#'                   }
#'                   }
#' @return The function returns a list containing the time series data generated. The list can contain multiple replicates (named as \code{Rep1}, \code{Rep2} etc.) equal to the \code{numReplicates} function argument.
#'         Each replicate can contain multiple targets (named as \code{Target1}, \code{Target2} etc.) based on the specified exposure space (\code{expSpace}). The \code{expSpace} and \code{controlFile} are also returned as part of this output list.
#' @seealso \code{createExpSpace}, \code{writeControlFile}, \code{viewModels}
#' @examples
#' # Example 1: Simple scaling
#' #-----------------------------------------------------------------------
#' attPerturb<-c("P_ann_tot_m","Temp_ann_avg_m")
#' attPerturbType = "regGrid"
#' attPerturbSamp = c(2, 2)
#' attPerturbMin = c(0.8, -1)
#' attPerturbMax = c(1.1, 1)
#' expSpace <- createExpSpace(attPerturb = attPerturb,
#'                            attPerturbSamp = attPerturbSamp,
#'                            attPerturbMin = attPerturbMin,
#'                            attPerturbMax = attPerturbMax,
#'                            attPerturbType = attPerturbType)
#' data(tankDat)
#' simScaling <- generateScenarios(reference = tank_obs,
#'                                 expSpace = expSpace,
#'                                 controlFile = "scaling")
#'
#' # Example 2: Seasonal scaling
#' #-----------------------------------------------------------------------
#' attPerturb<-c("P_ann_tot_m","P_ann_seasRatio")
#' attPerturbType = "regGrid"
#' attPerturbSamp = c(2, 2)
#' attPerturbMin = c(0.8, 0.9)
#' attPerturbMax = c(1.1, 1.2)
#' expSpace <- createExpSpace(attPerturb = attPerturb,
#'                            attPerturbSamp = attPerturbSamp,
#'                            attPerturbMin = attPerturbMin,
#'                            attPerturbMax = attPerturbMax,
#'                            attPerturbType = attPerturbType)
#' data(tankDat)
#' seasScaling <- generateScenarios(reference = tank_obs,
#'                                 expSpace = expSpace,
#'                                 controlFile = "scaling")
#'
#' # Example 3: Stochastic simulation using foreSIGHT default settings
#' #----------------------------------------------------------------------
#' \dontrun{
#' # create an exposure space
#' attPerturb <- c("P_ann_tot_m", "P_ann_nWet_m", "P_ann_R10_m")
#' attHold <- c("P_Feb_tot_m", "P_SON_dyWet_m", "P_JJA_avgWSD_m", "P_MAM_tot_m",
#' "P_DJF_avgDSD_m", "Temp_ann_rng_m", "Temp_ann_avg_m")
#' attPerturbType = "regGrid"
#' attPerturbSamp = c(2, 1, 1)
#' attPerturbMin = c(0.8, 1, 1)
#' attPerturbMax = c(1.1, 1, 1)
#' expSpace <- createExpSpace(attPerturb = attPerturb,
#'                            attPerturbSamp = attPerturbSamp,
#'                            attPerturbMin = attPerturbMin,
#'                            attPerturbMax = attPerturbMax,
#'                            attPerturbType = attPerturbType,
#'                            attHold = attHold,
#'                            attTargetsFile = NULL)
#' # load example data available in foreSIGHT
#' data(tankDat)
#' # perform stochastic simulation
#' simStochastic <- generateScenarios(reference = tank_obs,
#'                                    expSpace = expSpace,
#'                                    simLengthNyrs = 30)
#'                                    }
#' # Example 4: Simple Scaling with multi-site data
#' #-----------------------------------------------------------------------
#' attPerturb <- c("P_ann_tot_m","P_ann_seasRatio")
#' attPerturbType = "regGrid"
#' attPerturbSamp = c(3, 3)
#' attPerturbMin = c(0.8, 1.2)
#' attPerturbMax = c(0.8, 1.2)
#' expSpace <- createExpSpace(attPerturb = attPerturb,
#'                            attPerturbSamp = attPerturbSamp,
#'                            attPerturbMin = attPerturbMin,
#'                            attPerturbMax = attPerturbMax,
#'                            attPerturbType = attPerturbType)
#' # load multi-site rainfall data
#' data(barossaDat)
#' # perform simple scaling
#' simScaling <- generateScenarios(reference = barossa_obs,
#'                                 expSpace = expSpace,
#'                                 controlFile = "scaling")
#'
#' # Example 5: Multi-site stochastic simulation
#' #-----------------------------------------------------------------------
#' \dontrun{
#' attPerturb <- c("P_ann_tot_m")
#' attHold <- c("P_ann_wettest6monSeasRatio","P_ann_wettest6monPeakDay",
#' "P_ann_P99","P_ann_avgWSD_m", "P_ann_nWetT0.999_m")
#' attPerturbType = "regGrid"
#' # consider unperturbed climates in this example
#' attPerturbSamp = attPerturbMin = attPerturbMax = c(1)
#' expSpace <- createExpSpace(attPerturb = attPerturb,
#'                            attPerturbSamp = attPerturbSamp,
#'                            attPerturbMin = attPerturbMin,
#'                            attPerturbMax = attPerturbMax,
#'                            attPerturbType = attPerturbType,
#'                            attHold = attHold)
#' # load multi-site rainfall data
#' data(barossaDat)
#' # specify the penalty settings in a list
#' controlFileList <- list()
#' controlFileList[["penaltyAttributes"]] <- c("P_ann_tot_m",
#' "P_ann_wettest6monSeasRatio","P_ann_wettest6monPeakDay")
#' controlFileList[["penaltyWeights"]] <- c(0.5,0.5,0.5)
#' # specify the alternate model selections
#' controlFileList[["modelType"]] <- list()
#' controlFileList[["modelType"]][["P"]] <- "latent"
#' # specify model parameter selection
#' controlFileList[["modelParameterVariation"]] <- list()
#' controlFileList[["modelParameterVariation"]][["P"]] <- "harmonic"
#' # specify settings for multi-site model
#' controlFileList[["spatialOptions"]] <- list()
#' # specify spatial correlation perturbation factor
#' controlFileList[["spatialOptions"]][["spatCorFac"]] = 0.9
#' # write control file sttings to file
#' controlFileJSON <- jsonlite::toJSON(controlFileList, pretty = TRUE, auto_unbox = TRUE)
#' write(controlFileJSON, file = paste0(tempdir(), "controlFile.json"))
#' # run multi-site stochastic simulation - this will take a long time (e.g. hours)
#' sim <- generateScenarios(reference = barossa_obs, expSpace = expSpace,
#'                          controlFile = paste0(tempdir(), "controlFile.json"),seed=1)}
#' @export

generateScenarios <- function(reference,                # data frame of observed data with column names compulsary [$year, $month, $day, $P,] additional [$Temp, $RH, $PET, $uz, $Rs] (or a subset of these)
                              expSpace,           # the space contains multiple targets
                              simLengthNyrs = NULL,      # desired length of simulation in years
                              numReplicates = 1,  # reps
                              seedID = NULL,      # seed - user may set this to reproduce a previous simulation
                              controlFile = NULL   # NULL = default stochastic model options, "scaling" = simple scaling, json file = stochastic model options
                              ) {

  # Number of targets
  nTarget <- dim(expSpace$targetMat)[1]

  # Replicates and seed don't go with scaling
  if (!is.null(controlFile)) {
    if (controlFile == "scaling") {
      if (numReplicates > 1) stop("Simple scaling cannot generate replicates. Please set numReplicates to 1.")
      if (!is.null(seedID)) stop("Simple scaling cannot use a seed. Please set seedID to NULL.")
    }
  }

  # Create random seedID
  if (is.null(seedID)) {
    seedID <- round(stats::runif(1)*10000)
  }

  # Create seedID vector for all replicates
  if (numReplicates>0 & numReplicates%%1==0) {
    seedIDs <- seedID + seq(0, numReplicates-1)
    nRep <- length(seedIDs)
  } else {
    stop("numReplicates should be a positive integer")
  }

  allSim <- replicate(nRep, vector("list", nTarget), simplify = FALSE)

  for (iRep in 1:nRep) {

    cat(paste0("Generating replicate number ", iRep,  " out of ", nRep, " replicates...\n"))
    pb <- progress::progress_bar$new(
      #format = " [:bar] :elapsedfull",
      total = nTarget, clear = FALSE, width= 60)
    pb$tick(0)


    for (iTarg in 1:nTarget) {

      # Get the target location in the exposure space
      expTarg <- expSpace
      expTarg$targetMat <- expSpace$targetMat[iTarg, ]
      if(!is.null(expSpace$attRot)) {
        expTarg$attRot <- expSpace$attRot[iTarg]
      }


      #(    working on Target No. ", iTarg, " of ", nTarget, "\n"))
                 #"\n=============================================================\n"))
      pb$tick()
      # Call generateScenario for the target
      allSim[[iRep]][[iTarg]] <- generateScenario(reference = reference,
                                                  expTarg = expTarg,
                                                  simLengthNyrs = simLengthNyrs,
                                                  seedID = seedIDs[iRep],
                                                  controlFile = controlFile
                                                  )

      # Get & remove simDates and nml from the target simulation, will be added back later
      nmlOut <- allSim[[iRep]][[iTarg]][["nml"]]
      simDates <- allSim[[iRep]][[iTarg]][["simDates"]]
      allSim[[iRep]][[iTarg]][["nml"]] <- NULL
      allSim[[iRep]][[iTarg]][["simDates"]] <- NULL
    }
    names(allSim[[iRep]]) <- paste0("Target", 1:nTarget)
  }
  cat("Simulation completed")
  names(allSim) <- paste0("Rep", 1:nRep)
  allSim[["simDates"]] <- simDates
  allSim[["expSpace"]] <- expSpace
  allSim[["controlFile"]] <- nmlOut

  if (!is.null(controlFile)) {
    if (controlFile == "scaling") allSim[["controlFile"]] <- controlFile
  }

  return(allSim)

}


checkObsVars <- function(obs, file) {
  obsVars <- names(obs)[-which(names(obs) %in% c("year", "month", "day"))]
  for (o in obsVars) {
    if(!(o %in% fSVars)){
      logfile(paste0("Input variable ", o, " unrecognized."), file)
      logfile("Program terminated",file)
      stop(paste0("Input variable ", o, " unrecognized."))
    }
  }
  return(invisible())
}

# Process namelist
# reads controlFile and converts to required inputs: modelTag, modelInfoMon, attPenalty, optimArgs
# Anjana: pass "file" to readNamelist and modify checkNamelist functions to write errors into the file
getUserModelChoices <- function(controlFile, obs, attSel, file = NULL) {

  obsVars <- names(obs)[-which(names(obs) %in% c("year", "month", "day"))]
  attVars <- vapply(attSel,FUN = get.attribute.varType,FUN.VALUE=character(1),USE.NAMES = FALSE)

  if (!is.null(file)) {
    checkObsVars(obs, file)

    for (i in 1:length(attSel)) {
        if(!(attVars[i] %in% obsVars)){
          logfile(paste0("Observations do not contain the variable ", attVars[i], " to compute the selected attribute attSel [",i,"]"), file)
          logfile("Program terminated",file)
          stop(paste0("Observations do not contain the variable ", attVars[i], " to compute the selected attribute attSel [",i,"]"))
        }
    }
  }

  # defaults
  modelInfoMod <- list()
  attPenalty <- NULL
  optimArgs <- list()
  spatialArgs = list()

  if (is.null(controlFile)) {
    modelTag <- NULL
    for (v in intersect(obsVars, attVars)) {
      modelTag <- c(modelTag, getModelTag(nml= NULL, v))
    }
  } else if (controlFile == "scaling") {
      is.seas = grepl('seasRatio',attSel)
      len.seas = length(which(is.seas))
      if (len.seas==0){
        modelTag = c('Simple-ann')
      } else if (length(attSel)==2*len.seas){
        modelTag = c('Simple-seas')
      } else {
        modelTag = c('Simple-ann','Simple-seas')
      }
    } else {

    # Read in user namelist
    # Some observations may have a user choice - others may be default
    # It is assumed that if the input obs contains a variable it is to be simulated provided that a correspondiong attribute is selected

    # JSON file, else it must alread by an Rlist
    if (is.character(controlFile)) {
      nml <- readNamelist(controlFile)
    } else if (is.list(controlFile)) {
      nml <- controlFile
    }
    nmlVars <- names(nml[["modelType"]])
    allVars <- union(unique(attVars), nmlVars)

    modelTag <- NULL
    for (v in allVars) {
      modelTag <- c(modelTag, getModelTag(nml = nml, v))
    }
    modelInfoMod <- getModelInfoMod(nml)
    optimArgs <- getOptimArgs(nml)
    attPenalty <- getAttPenalty(nml)
    spatialArgs = list(sites=nml$spatialOptions$sites,spatCorMatIn=nml$spatialOptions$spatCorMatIn,spatCorFac=nml$spatialOptions$spatCorFac)
  }

  return(list(modelTag = modelTag,
              modelInfoMod = modelInfoMod,
              optimArgs = optimArgs,
              attPenalty = attPenalty,
              spatialArgs = spatialArgs))

}

# add extra model info for simple/seasonal scaling (modelinfo[[mod]]$simvar), and perform some extra checks
# note: this could possibly be incorporated into get.multi.model.info or check_duplicates_mismatch but would significant changes to either function
add_scaling_info = function(obs,attSel,modelInfo){

  obsVars <- names(obs)[-which(names(obs) %in% c("year", "month", "day"))]
  attVars <- vapply(attSel,FUN = get.attribute.varType,FUN.VALUE=character(1),USE.NAMES = FALSE)
  for (v in intersect(obsVars, attVars)) {
    attSelVar = attSel[attVars==v]
    tmp = unlist(strsplit(attSelVar,paste0(v,'_')))
    suffix = tmp[tmp!='']
    if (any(suffix%in%c('ann_tot_m','ann_avg_m'))){
      if (length(attSelVar)==1){
        modelInfo[["Simple-ann"]]$simVar = c(modelInfo[["Simple-ann"]]$simVar,v)
      } else if ((length(attSelVar)==2)){
        if (any(grepl('seasRatio',attSelVar))){
          modelInfo[["Simple-seas"]]$simVar = c(modelInfo[["Simple-seas"]]$simVar,v)
        }
      } else {
        stop("Scaling only works with 2 attributes (total/avg and seasRatio) for ",v)
      }
    } else {
      stop("Require attribute with name '",v,"_ann_tot_m' or ",v,"_ann_avg_m' for scaling")
    }
  }

  return(modelInfo)
}

#' Produces time series of hydroclimatic variables for an exposure target.
#'
#' \code{generateScenario} is the base function used by \code{generateScenarios}.
#' The function produces time series of hydroclimatic variables using requested climate attributes that correspond to a single target in the exposure space.
#' The function argument definitions are detailed in the documentation of \code{generateScenarios}; please refer to that documentation using \code{?generateScenarios}.
#' @inheritParams generateScenarios
#' @param expTarg a named vector; the attributes at the target location in the exposure space
#' \code{generateScenario} is intended to be used to adapt the functionality of \code{generateScenarios} for use in a parallel computing environment.
#' @seealso \code{generateScenarios}
#' @export
#' @import GA

generateScenario <- function(reference,       # data frame of observed data with column names compulsary [$year, $month, $day, $P,] additional [$Temp, $RH, $PET, $uz, $Rs] (or a subset of these)
                             expTarg,
                             simLengthNyrs = NULL,
                             seedID = NULL,
                             controlFile = NULL
){

  # renamed obs to reference (rename everywhere sometime)
  obs <- reference

  # Create random seedID
  if (is.null(seedID)) {
    seedID <- round(stats::runif(1)*10000)
  }
  set.seed(seedID)

  file <- paste0(tempdir(), "/generateScenario_log.txt")

  # Checking
  #-------------------------------------------------------
  banner("CHECK FOR ARGUMENT INPUTS",file)
  progress("Checking argument inputs...",file)


  # Unpack expTarg
  #------------------------------------------------------
  attPerturb <- expTarg$attPerturb
  attHold <- expTarg$attHold
  attSel <- c(attPerturb,attHold)

  # Process Namelist
  #------------------------------------------------------
  userModelChoices <- getUserModelChoices(controlFile, obs, attSel, file)
  modelTag <- userModelChoices$modelTag
  modelInfoMod <- userModelChoices$modelInfoMod
  optimArgs <- userModelChoices$optimArgs
  attPrim <- userModelChoices$attPenalty
  spatialArgs = userModelChoices$spatialArgs

  # Update optimArgs
  #------------------------------------------------------
  optimArgs <- utils::modifyList(optimArgsdefault, optimArgs)


  # Identify target
  #-------------------------------------------------------
  if(!is.null(dim(expTarg$targetMat))) {
    if (dim(expTarg$targetMat)[1] > 1) {
      progress("expTarg has multiple targets. 'generateScenario' works on a single target, use 'generateScenarios' instead", file)
      stop()
    }
  }
  targetMat <- expTarg$targetMat

  check_duplicates_mismatch(obs=obs,
                            attSel=attSel,
                            attPrim=attPrim,
                            attHold=attHold,
                            attPerturb=attPerturb,
                            modelTag=modelTag,
                            optimArgs=optimArgs,
                            file=file)
  progress("Argument input format OK",file)

  banner("CHECK FOR MODEL AND ATTRIBUTE COMBINATIONS",file)
  progress("Checking model and attribute combinations...",file)

  check_models_attributes(names=names(obs),
                          attSel=attSel,
                          attPrim=attPrim,
                          modelTag=modelTag,
                          file=file)

  progress("Model and attribute combinations OK",file)

  #CHECK FOR INPUTS
  banner("CHECK FOR DATAFRAME INPUT",file)
  progress("Checking dataframe input...",file)
  inputcheck<-input_check(obs,file,simLengthNyrs)
  obs=inputcheck$data                                      # USE NEW APPENDED/CHECKED DATA
  progress("Dataframe input OK",file)

  #GET ADDITIONAL MODEL INFO, ATT INFO & SORT (make into separate script/functions)
  nMod=length(modelTag)
  modelInfo=get.multi.model.info(modelTag=modelTag)

  # add extra model info for simple/seasonal scaling and perform some extra checks
  if(modelTag[1]%in%c("Simple-ann","Simple-seas")){
    modelInfo = add_scaling_info(obs=obs,attSel=attSel,modelInfo=modelInfo)
  }

  #UPDATE MODELINFO IF NEEDED
  for(mod in 1:nMod){
    if(!is.null(modelInfoMod[[modelTag[mod]]])){
      #modifyList
      if(mod==1) progress("Updating model info...",file)
      defaultMods=list(minBound=NULL,maxBound=NULL,fixedPars=NULL)
      modPars=utils::modifyList(defaultMods,modelInfoMod[[modelTag[mod]]])
      modelInfo[[modelTag[mod]]]=update.model.info(modelTag=modelTag[mod],
                                                   modelInfo=modelInfo[[modelTag[mod]]],
                                                   fixedPars=modPars$fixedPars,
                                                   minUserBound=modPars$minBound,
                                                   maxUserBound=modPars$maxBound,
                                                   file=file)  #need to build in checks for this
    }
  }

  modelTag=update.simPriority(modelInfo=modelInfo)
  simVar=sapply(X=modelInfo[modelTag],FUN=return.simVar,USE.NAMES=TRUE)       #?CREATE MODEL MASTER INFO - HIGHER LEVEL?

  attInfo=attribute.info.check(attSel=attSel,attPrim=attPrim, lambda.mult = optimArgs$lambda.mult)                                 # vector of selected attributes (strings)
  if(modelTag[1]%in%c("Simple-ann","Simple-seas")){simVar=attInfo$varType}
  attInd=get.att.ind(attInfo=attInfo,simVar=simVar)
  attInfo=update.att.Info(attInfo=attInfo,attInd=attInd,modelTag=modelTag,simVar=simVar) #add extra level for easier model mangmt
  if(!(modelTag[1]%in%c("Simple-ann","Simple-seas"))){
    nParTot=0
    for(i in 1:nMod){
      nParTot=nParTot+modelInfo[[i]]$npars #total number of pars
    }
  }

  #GET DATES DATA (and indexes for harmonic periods)
  banner("INDEXING DATES",file)
  progress("Indexing dates...",file)
  dateExtnd=dateExtender(obs=obs,simLengthNyrs=simLengthNyrs,file=file,modelTag=modelTag)  #Extend dates if needed
  datInd=mod.get.date.ind.extnd(obs=obs,dateExtnd=dateExtnd,modelTag=modelTag,modelInfo=modelInfo,simLengthNyrs=simLengthNyrs,file=file,southHemi=TRUE)

  progress("Dates indexed OK",file)

  # Adding AR1 parameter for P-har-wgen to ModelInfo (Anjana: moved to a separate fn from modelSequencer)
  # Has to be after datInd since dates are used in the AR1 param calculation
  modelInfo <- add_ar1Param(modelTag = modelTag,
                            modelInfo = modelInfo,
                            datInd = datInd)

  #-------------------------------------
  #SIMPLE/SEASONAL SCALING
  if(modelTag[[1]]%in%c("Simple-ann","Simple-seas")){
    sim = list()
    for (mod in 1:nMod){
      if (modelTag[mod]=='Simple-ann'){
        # simple scaling
        banner("SIMPLE SCALING FOR OBSERVED DATA",file)
        progress("Simple scaling data...",file)
        i_simple_ann = which(attInfo$varType %in% modelInfo[["Simple-ann"]]$simVar)
        varSel = modelInfo[["Simple-ann"]]$simVar[i_simple_ann]
        sim[varSel]=simple.scaling(target=unlist(targetMat)[i_simple_ann],
                                   targetType=attInfo$targetType[i_simple_ann],
                                   data=obs,
                                   varType=attInfo$varType[i_simple_ann],
                                   period=modelInfo[[modelTag[mod]]]$nperiod,
                                   i.pp=datInd[[modelTag[mod]]]$i.pp)[varSel]
        progress("Simple scaling OK",file)
      } else if (modelTag[mod]=='Simple-seas'){
        # seasonal scaling
        banner("SEASONAL SCALING FOR OBSERVED DATA",file)
        progress("Seasonal scaling data...",file)
        for (v in modelInfo[["Simple-seas"]]$simVar){
          i1 = which(names(targetMat)%in%paste0(v,c('_ann_tot_m','_ann_avg_m'))) # tot/avg attribute
          i2 = which(grepl(paste0(v,'_ann_seasRatio'),names(targetMat)))         # seasonality ratio attribute
          attSeas = names(targetMat)[i2]
          o = attribute.calculator.setup(attSeas,datInd[[modelTag[mod]]])
          sim[[v]] = simple.scaling.seasonal(target_total_fac=unlist(targetMat)[i1],
                                             target_seas_fac=unlist(targetMat)[i2],
                                             data=obs,
                                             varType=v,
                                             i.T= datInd[[modelTag[mod]]]$i.pp[[1]],
                                             i.S2 = o[[attSeas]]$attArgs$indexWet,
                                             i.S1=o[[attSeas]]$attArgs$indexDry,
                                             phi=o[[attSeas]]$attArgs$phi)[[v]]
        }
        progress("Seasonal scaling OK",file)
      }
    }
    sim[["simDates"]] <- dateExtnd
    progress("Data scaled OK",file)

  }else{

    #STOCHASTIC SIMULATION
    #---------------------------

    nsite = dim(obs[[simVar[1]]])[2] ## need to fix this up for nMod > 1

    if (nsite==1){

      #GET ATTRIBUTES OF OBSERVED DATA (testing attribute calc function)
      banner("OBSERVED BASELINE ATTRIBUTE CALCULATION",file)
      progress("Calculating attributes...",file)

      attObs=list()                            #make this into its own function (also inserted into model sequencer)
      for(i in 1:nMod){
        attObs[[i]]=attribute.calculator(attSel=attSel[attInd[[i]]],
                                         data=obs[[simVar[i]]],
                                         datInd=datInd[["obs"]])#,
        #                                             attribute.funcs=attribute.funcs)
      }
      attObs=unlist(attObs); attObs=attObs[attSel]   #unlist attObs and make sure order is correct

      progress(paste("Attributes of observed series - ",paste(attSel,": ",signif(attObs,digits=5),collapse = ", ",sep=""),sep=""),file)
      progress("Attributes calculated OK",file)   #NEED SOME ACTUAL CHECKING HERE BEFORE PRONOUNCING OK

      #OPTIMISING TO DETERMINE PARS
      #LOOP OVER EXPOSURE SPACE POINTS TO DETERMINE PARS
      banner("DETERMINING STOCHASTIC MODEL PARAMETERS FOR TARGET",file)
      progress("Determining stochastic model parameters at target location...",file)
      progress("Simulating stochastic time series...",file)
      progress("Starting cluster...",file)

      #DETERMINE WHICH PARS ATTACH TO WHICH MODEL (make this a function in stochParManager.R)
      parLoc=whichPars(simVar=simVar,modelInfo=modelInfo)

      #SCREEN INAPPROPRIATE SUGGESTIONS IF ANY
      if(!is.null(optimArgs$suggestions)){
        optimArgs$suggestions=screenSuggest(suggest=optimArgs$suggestions,modelInfo=modelInfo,modelTag=modelTag,parLoc=parLoc)
      }

      a<-Sys.time()

      #IF "OAT" ROTATE attPrim
      #attRot is returned only for "OAT" grids
      if(!is.null(expTarg$attRot)) {
        attApp <- expTarg$attRot
      } else {
        attApp <- attPrim
      }

      sim=simulateTarget(optimArgs=optimArgs,         #sim[[i]]$P, $Temp $attSim $targetSim
                         simVar=simVar,
                         modelTag=modelTag,
                         modelInfo=modelInfo,
                         attSel=attSel,
                         attPrim=attApp,              #controlled via switch
                         attInfo=attInfo,
                         attInd=attInd,
                         datInd=datInd,
                         initCalibPars=NULL,
                         targetLoc=targetMat,     #is  a vector  (just 1 target here)
                         attObs=attObs,
                         parLoc=parLoc,
                         parSim=NULL,
                         # Anjana - do I need this?
                         setSeed=seedID,                   #seed based on loop counter
                         file=file)

      b<-Sys.time()
      runt <- b-a
      logfile(signif(runt,digits=3),file)
      nmlOut <- toNamelist(modelTag = modelTag, modelInfoMod = modelInfoMod, optimArgs = optimArgs, attPenalty = attPrim)
      sim[["nml"]] <- nmlOut
      sim[["simDates"]] <- dateExtnd
      progress("Stochastic model parameters and time series obtained at target location ", file)

    } else {

      #OPTIMISING TO DETERMINE PARS
      #LOOP OVER EXPOSURE SPACE POINTS TO DETERMINE PARS
      banner("DETERMINING STOCHASTIC MODEL PARAMETERS FOR TARGET",file)
      progress("Determining stochastic model parameters at target location...",file)
      progress("Simulating stochastic time series...",file)
      progress("Starting cluster...",file)

      #DETERMINE WHICH PARS ATTACH TO WHICH MODEL (make this a function in stochParManager.R)
      parLoc=whichPars(simVar=simVar,modelInfo=modelInfo)

      #SCREEN INAPPROPRIATE SUGGESTIONS IF ANY
      if(!is.null(optimArgs$suggestions)){
        optimArgs$suggestions=screenSuggest(suggest=optimArgs$suggestions,modelInfo=modelInfo,modelTag=modelTag,parLoc=parLoc)
      }

      a<-Sys.time()

      #IF "OAT" ROTATE attPrim
      #attRot is returned only for "OAT" grids
      if(!is.null(expTarg$attRot)) {
        attApp <- expTarg$attRot
      } else {
        attApp <- attPrim
      }

      # perform Stage 1 of multi-site simulation, where we calculate marginal parameters at each site, based on input spatial correlation matrix (spatialArgs$spatCorMatIn)
      sim1=simulateTargetMarg(optimArgs=optimArgs,
                              simVar=simVar,
                              modelTag=modelTag,
                              modelInfo=modelInfo,
                              attSel=attSel,
                              attPrim=attApp,
                              attInfo=attInfo,
                              attInd=attInd,
                              datInd=datInd,
                              initCalibPars=NULL,
                              targetLoc=targetMat,
                              parLoc=parLoc,
                              parSim=NULL,
                              setSeed=seedID,
                              file=file,
                              nMod=nMod,
                              obs=obs,
                              spatialArgs=spatialArgs)

      # perform Stage 2 of multi-site simulation, using marginal parameters from Stage 1 and this time calculating spatial correlation matrix
      sim2 = simulateTargetCor(optimArgs=optimArgs,
                               simIn=sim1,
                               simVar=simVar,
                               modelTag=modelTag,
                               modelInfo=modelInfo,
                               attSel=attSel,
                               attPrim=attApp,
                               attInfo=attInfo,
                               attInd=attInd,
                               datInd=datInd,
                               targetLoc=targetMat,
                               parLoc=parLoc,
                               setSeed=seedID,
                               file=file,
                               nMod=nMod,
                               obs=obs,
                               spatialArgs=spatialArgs)

      # perform Stage 3 of multi-site simulation, where we re-calculate marginal parameters at each site, based on spatial correlation matrix from Stage 2
      spatialArgs3 = spatialArgs
      spatialArgs3$spatCorMatIn = sim2$cor_par
      sim3=simulateTargetMarg(optimArgs=optimArgs,
                              simVar=simVar,
                              modelTag=modelTag,
                              modelInfo=modelInfo,
                              attSel=attSel,
                              attPrim=attApp,
                              attInfo=attInfo,
                              attInd=attInd,
                              datInd=datInd,
                              initCalibPars=NULL,
                              targetLoc=targetMat,
                              parLoc=parLoc,
                              parSim=NULL,
                              setSeed=seedID,
                              file=file,
                              nMod=nMod,
                              obs=obs,
                              spatialArgs=spatialArgs3)

      sim = list()
      sim$Stage1 = sim1
      sim$Stage2 = sim2
      sim$Stage3 = sim3

      b<-Sys.time()
      runt <- b-a
      logfile(signif(runt,digits=3),file)
      nmlOut <- toNamelist(modelTag = modelTag, modelInfoMod = modelInfoMod, optimArgs = optimArgs, attPenalty = attPrim)
      sim[["nml"]] <- nmlOut
      sim[["simDates"]] <- dateExtnd
      progress("Stochastic model parameters and time series obtained at target location ", file)

    }

  } #END STOCHASTIC SEGMENT

  return(sim)
}

#----------------------------------------------------------------------------

simulateTargetMarg = function(optimArgs=NULL,
                              simVar=NULL,
                              modelTag=NULL,
                              modelInfo=NULL,
                              attSel=NULL,
                              attPrim=NULL,
                              attInfo=NULL,
                              attInd=NULL,
                              datInd=NULL,
                              initCalibPars=NULL,
                              targetLoc=NULL,
                              parLoc=NULL,
                              parSim=NULL,
                              setSeed=NULL,
                              file=NULL,
                              nMod=NULL,
                              obs=NULL,
                              spatialArgs=NULL){

  nday = datInd$obs$ndays

  ### DM NOTE: CODE CURRENTLY NOT PROPERLY SETUP TO WORK WITH MULTIPLE MODELS - PROB NEED TO ADD MODEL TO SIM LIST
  sim = list(sites=NULL)

  for(i in 1:nMod){
    nsite = dim(obs[[simVar[i]]])[2]
    sites = colnames(obs[[simVar[i]]])
    if (is.null(sites)){sites=paste0('site',seq(1,nsite))}
    if ((nsite>1)&(nMod>1)){
      stop('cannot have (nsite>1)&(nMod>1)')
    }

    set.seed(setSeed)
    if (nsite==1){
      spatCorMatIn = NULL
      MVTsampleMat = matrix(stats::rnorm(n=nday,mean=0.,sd=1.),ncol=1)
    } else {
      if (!is.null(spatialArgs$spatCorMatIn)){
        spatCorMatIn = spatialArgs$spatCorMatIn
      } else {
        spatCorMatIn = diag(nsite)
      }
#      MVTsampleMat = mvtnorm::rmvnorm(n=nday,sigma=spatCorMatIn)
    spatCorMatIn_PD = as.matrix(Matrix::nearPD(spatCorMatIn,keepDiag = T)$mat)
    MVTsampleMat = mvtnorm::rmvnorm(n=nday,sigma=spatCorMatIn_PD)

      colnames(MVTsampleMat) = sites
    }

    simMultiSite = list(sites=NULL,P=NULL)  ### NEED TO FIX THIS UP FOR VARIABLES OTHER THAN P
    for (s in 1:nsite){
      site = sites[s]
      cat(site,'\n')
      #GET ATTRIBUTES OF OBSERVED DATA (testing attribute calc function)
      banner("OBSERVED BASELINE ATTRIBUTE CALCULATION",file)
      progress("Calculating attributes...",file)
      attObs=attribute.calculator(attSel=attSel[attInd[[i]]],
                                     data=obs[[simVar[i]]][,s],
                                     datInd=datInd[["obs"]])

      attObs=unlist(attObs); attObs=attObs[attSel]   #unlist attObs and make sure order is correct

      progress(paste("Attributes of observed series - ",paste(attSel,": ",signif(attObs,digits=5),collapse = ", ",sep=""),sep=""),file)
      progress("Attributes calculated OK",file)   #NEED SOME ACTUAL CHECKING HERE BEFORE PRONOUNCING OK

      simMultiSite$sites[[site]] = simulateTarget(optimArgs=optimArgs,
                                         simVar=simVar,
                                         modelTag=modelTag,
                                         modelInfo=modelInfo,
                                         attSel=attSel,
                                         attPrim=attPrim,
                                         attInfo=attInfo,
                                         attInd=attInd,
                                         datInd=datInd,
                                         initCalibPars=initCalibPars,
                                         targetLoc=targetLoc,
                                         attObs=attObs,
                                         parLoc=parLoc,
                                         parSim=parSim,
                                         setSeed=setSeed,
                                         file=file,
                                         randomUnitNormalVector=MVTsampleMat[,s])

      if (nsite>1){simMultiSite$P$sim = cbind(simMultiSite$P$sim,simMultiSite$sites[[site]]$P$sim)}

    }
  }

  if (nsite>1){
    sim = simMultiSite
    sim$cor_par = spatCorMatIn
  } else {
    sim = simMultiSite$sites[[1]]
  }
  return(sim)

}

#----------------------------------------------------------------------------

simulateTargetCor = function(optimArgs=NULL,
                             simIn=NULL,
                             simVar=NULL,
                             modelTag=NULL,
                             modelInfo=NULL,
                             attSel=NULL,
                             attPrim=NULL,
                             attInfo=NULL,
                             attInd=NULL,
                             datInd=NULL,
                             targetLoc=NULL,
                             attObs=NULL,
                             parLoc=NULL,
                             setSeed=NULL,
                             file=NULL,
                             randomUnitNormalVector=NULL,  ### remove this
                             nMod=NULL,
                             obs=NULL,
                             spatialArgs=NULL){

  if (nMod>1){
    stop('cannot have nMod>1')
  } else {
    mod = 1
  }
  sites = names(simIn$sites)
  nsite=length(sites)
  nday = datInd$obs$ndays

  modelInfoSites = list()
  for (s in 1:nsite){
    site = sites[s]
    modelInfoSites[[site]] = modelInfo
    modelInfoSites[[site]][[modelTag[mod]]]$minBound = modelInfoSites[[site]][[modelTag[mod]]]$maxBound = simIn$sites[[site]]$parS
  }

  cor_par = matrix(nrow=nsite,ncol=nsite)
  for (s1 in 1:(nsite-1)){
    site1 = sites[s1]
    cor_par[s1,s1] = 1

    for (s2 in (s1+1):nsite){
      site2 = sites[s2]
      cor_obs = stats::cor(obs$P[,c(s1,s2)])[1,2]
      if (!is.null(spatialArgs$spatCorFac)){cor_obs=spatialArgs$spatCorFac*cor_obs}

      rho_1_2_list = seq(0.01,1,0.01)
      cor_sim_list = c()

      corMat = matrix(nrow=2,ncol=2)
      corMat[1,1] = corMat[2,2] = 1

      for (i in 1:length(rho_1_2_list)){

        corMat[1,2] = corMat[2,1] = rho_1_2_list[i]

        set.seed(setSeed)
#        MVTsampleMat = mvtnorm::rmvnorm(n=nday,sigma=corMat)
        corMat_PD = as.matrix(Matrix::nearPD(corMat,keepDiag = T)$mat)
        MVTsampleMat = mvtnorm::rmvnorm(n=nday,sigma=corMat_PD)

        attObs1=attribute.calculator(attSel=attSel[attInd[[mod]]],
                                    data=obs[[simVar[mod]]][,s1],
                                    datInd=datInd[["obs"]])
        attObs1=unlist(attObs1); attObs1=attObs1[attSel]
        sim1 = simulateTarget(optimArgs=optimArgs,
                              simVar=simVar,
                              modelTag=modelTag,
                              modelInfo=modelInfoSites[[site1]],
                              attSel=attSel,
                              attPrim=attPrim,
                              attInfo=attInfo,
                              attInd=attInd,
                              datInd=datInd,
                              targetLoc=targetLoc,
                              attObs=attObs1,
                              parLoc=parLoc,
                              file=file,
                              randomUnitNormalVector=MVTsampleMat[,1])

        attObs2=attribute.calculator(attSel=attSel[attInd[[mod]]],
                                     data=obs[[simVar[mod]]][,s2],
                                     datInd=datInd[["obs"]])
        attObs2=unlist(attObs2); attObs2=attObs2[attSel]
        sim2 = simulateTarget(optimArgs=optimArgs,
                              simVar=simVar,
                              modelTag=modelTag,
                              modelInfo=modelInfoSites[[site2]],
                              attSel=attSel,
                              attPrim=attPrim,
                              attInfo=attInfo,
                              attInd=attInd,
                              datInd=datInd,
                              targetLoc=targetLoc,
                              attObs=attObs2,
                              parLoc=parLoc,
                              file=file,
                              randomUnitNormalVector=MVTsampleMat[,2])

        cor_sim_list[i] = stats::cor(sim1$P$sim,sim2$P$sim)

      }

      diff = abs(cor_sim_list-cor_obs)
      cor_par[s1,s2] = cor_par[s2,s1] = rho_1_2_list[which(diff==min(diff))]
    }
    cor_par[nsite,nsite] = 1

  }

  sim = list(sites=NULL,P=NULL)
  set.seed(setSeed)
#  MVTsampleMat = mvtnorm::rmvnorm(n=nday,sigma=cor_par)
  cor_par_PD = as.matrix(Matrix::nearPD(cor_par,keepDiag = T)$mat)
  MVTsampleMat = mvtnorm::rmvnorm(n=nday,sigma=cor_par_PD)

  for (s in 1:nsite){
    site = sites[s]
    sim$sites[[site]] = simulateTarget(optimArgs=optimArgs,
                                 simVar=simVar,
                                 modelTag=modelTag,
                                 modelInfo=modelInfoSites[[site]],
                                 attSel=attSel,
                                 attPrim=attPrim,
                                 attInfo=attInfo,
                                 attInd=attInd,
                                 datInd=datInd,
                                 targetLoc=targetLoc,
                                 attObs=attObs2,
                                 parLoc=parLoc,
                                 file=file,
                                 randomUnitNormalVector=MVTsampleMat[,s])

    sim$P$sim = cbind(sim$P$sim,sim$sites[[site]]$P$sim)

  }

  sim$cor_par = cor_par

  return(sim)

}

#----------------------------------------------------------------------------

# writeScenarioOut_temp <- function(){
#
#
#   # #WRITE SIMPLE SCALING OUTPUT TO CSV & RDATA FILE
#   # if(IOmode=="verbose"){
#   #   if(!isTRUE(dir.exists(paths$Scenario))){
#   #     dir.create(paths$Scenario)
#   #   }
#   #   if(!isTRUE(dir.exists(paths$CSV))){
#   #     dir.create(paths$CSV)
#   #   }
#   #   if(!isTRUE(dir.exists(paths$RData))){
#   #     dir.create(paths$RData)
#   #   }
#   #
#   #   fnam<-simpleSaveTarget(data=sim,
#   #                     dates=obs[,c("year","month","day")],
#   #                     simVar=simVar,
#   #                     attSel=attSel,
#   #                     target=targetMat[i,],
#   #                     modelTag=modelTag[1],
#   #                     modelInfo=modelInfo[1],
#   #                     paths=paths)
#   # }
#
#   if(IOmode == "dev"){
#     #**DEV MODE TOWN**
#     #MATRICES BECOME VECTORS (1 TARGET AT A TIME)
#     simTarget=sim[[i]]$targetSim
#     parSim=sim[[i]]$parS
#     attValue=sim[[i]]$attSim  #NEW
#     req=targetMat[i,]
#     objScore=sim[[i]]$score
#
#     simTarget[length(attSel)+1]=parSim[nParTot+1]=attValue[length(attSel)+1]=req[length(attSel)+1]=arrayID
#     simTarget[length(attSel)+2]=parSim[nParTot+2]=attValue[length(attSel)+2]=req[length(attSel)+2]=seedCatalogue[i]
#
#     targ=t(matrix(simTarget))
#     value=t(matrix(attValue))
#     param=t(matrix(parSim))
#     request=t(matrix(req))
#     fitnesses=t(matrix(objScore))
#
#     #WRITE SIMS TO CSV FILE
#     write.table(targ,paste0("target",i,".csv"),append=FALSE,col.names=FALSE,row.names=FALSE,sep=",")
#     write.table(value,paste0("value",i,".csv"),append=FALSE,col.names=FALSE,row.names=FALSE,sep=",")
#     write.table(param,paste0("pars",i,".csv"),append=FALSE,col.names=FALSE,row.names=FALSE,sep=",")
#     write.table(request,paste0("request",i,".csv"),append=FALSE,col.names=FALSE,row.names=FALSE,sep=",")
#     write.table(fitnesses,paste0("fitness",i,".csv"),append=FALSE,col.names=FALSE,row.names=FALSE,sep=",")
#
#     simDat=makeOutputDataframe(data=sim[[i]],dates=dateExtnd,simVar=simVar,modelTag=modelTag[1])
#     write.table(simDat,file=paste0("scenario",i,".csv"),row.names=FALSE,quote = FALSE,sep=",")
#
#     devPlotSummary(obs=obs,
#                    sim=sim[[i]],
#                    simVar=simVar,
#                    datInd=datInd,
#                    attSel=attSel,
#                    attPrim=attPrim,
#                    simTarget=simTarget[c(1:length(attSel))],  #will be different for ... non arrayID
#                    target=targetMat[i,],
#                    targetType=attInfo$targetType,
#                    modelTag=modelTag,
#                    optimArgs=optimArgs,
#                    id=i,
#                    nTarget=nTarget,
#                    IOmode=IOmode)
#   }else if (IOmode=="verbose"){
#
#     if(!isTRUE(dir.exists(paths$Scenario))){
#       dir.create(paths$Scenario)
#     }
#     if(!isTRUE(dir.exists(paths$CSV))){
#       dir.create(paths$CSV)
#     }
#     if(!isTRUE(dir.exists(paths$RData))){
#       dir.create(paths$RData)
#     }
#     #SAVE EVERYTHING
#     simTarget[i,]=sim[[i]]$targetSim
#     parSim[i,]=sim[[i]]$parS
#     attValue[i,]=sim[[i]]$attSim  #NEW
#
#     #STORE IN RDATA &  WRITE DATA FRAME TO CSV
#     progress(p("Writing scenario ",i," to .RData and .csv"),file)
#     fnam[i]=saveTarget(data=sim[[i]],
#                        dates=dateExtnd,
#                        modelTag=modelTag,
#                        modelInfo=modelInfo,
#                        simVar=simVar,
#                        target=targetMat[i,],
#                        attSel=attSel,
#                        attPrim=attPrim,
#                        paths=paths)
#
#     if(!isTRUE(dir.exists(paths$Diagnostics))){
#       dir.create(paths$Diagnostics)
#     }#OUTPUT PLOT SUMMARY   #need alt filename here
#
#     devPlotSummary(obs=obs,
#                    sim=sim[[i]],
#                    simVar=simVar,
#                    datInd=datInd,
#                    attSel=attSel,
#                    attPrim=attPrim,
#                    simTarget=simTarget[i,c(1:length(attSel))],
#                    target=targetMat[i,],
#                    targetType=attInfo$targetType,
#                    modelTag=modelTag,
#                    optimArgs=optimArgs,
#                    id=i,
#                    nTarget=nTarget,
#                    IOmode=IOmode,
#                    paths=paths)
#
#   } else {
#     progress("No output generated in suppress mode",file)
#   }
#
# }


# writeScenarioOut <- function(sim = NULL,
#                              dates = NULL,
#                              simVar = NULL,
#                              attSel = NULL,
#                              target = NULL,
#                              modelTag = NULL,
#                              modelInfo = list(),
#                              outputPaths = NULL,
#                              suffixFileName = NULL
# ){
#
#   if (modelTag == "Simple-ann") {
#
#     fnam<-simpleSaveTarget(data=sim,
#                            dates=dates,
#                            simVar=simVar,
#                            attSel=attSel,
#                            target=targetMat,
#                            modelTag=modelTag,
#                            modelInfo=modelInfo,
#                            paths=paths,
#                            suffixFileName = suffixFileName)
#   } else {
#
#     simTarget=sim$targetSim
#     parSim=sim$parS
#     attValue=sim$attSim  #NEW
#     req=targetMat
#     objScore=sim$score
#
#
#     targ=t(matrix(simTarget))
#     value=t(matrix(attValue))
#     param=t(matrix(parSim))
#     request=t(matrix(req))
#     fitnesses=t(matrix(objScore))
#
#     #WRITE SIMS TO CSV FILE
#     write.table(targ,paste0("target",i,".csv"),append=FALSE,col.names=FALSE,row.names=FALSE,sep=",")
#     write.table(value,paste0("value",i,".csv"),append=FALSE,col.names=FALSE,row.names=FALSE,sep=",")
#     write.table(param,paste0("pars",i,".csv"),append=FALSE,col.names=FALSE,row.names=FALSE,sep=",")
#     write.table(request,paste0("request",i,".csv"),append=FALSE,col.names=FALSE,row.names=FALSE,sep=",")
#     write.table(fitnesses,paste0("fitness",i,".csv"),append=FALSE,col.names=FALSE,row.names=FALSE,sep=",")
#
#     simDat=makeOutputDataframe(data=sim[[i]],dates=dateExtnd,simVar=simVar,modelTag=modelTag[1])
#     write.table(simDat,file=paste0("sim_",suffixFileName,".csv"),row.names=FALSE,quote = FALSE,sep=",")
#
#   }
#
# }


# plotDiagnosticsx <- function() {
#
#   #---------------------------------------------------------
#   #SIMULATE AUXILLIARY TS
#   #PET and any other auxilliary vars
#   banner("SUMMARISING SIMULATION DIAGNOSTICS",file)
#   #--------------------------------------------------------
#   #SIMULATION INFORMATION SAVING
#    if(IOmode == "verbose"){  #IF RUN IN VERBOSE MODE
#      if(!isTRUE(dir.exists(paths$Metadata))){
#         dir.create(paths$Metadata)
#      }
#      progress("Saving simulation information to file - requested targets, simulated targets, optimised parameters",file)
#      if(modelTag[1] != "Simple-ann"){
#        write.table(simTarget,file=paste0(paths$Metadata,"/simulatedTargets.csv"),row.names=FALSE,quote = FALSE,sep=",",col.names=attSel)     #WRITE SIMULATED TARGETS TO CSV
#        parNames=NULL; for(i in 1:length(simVar)){parNames=c(parNames,modelInfo[[i]]$parNam)}
#        write.table(parSim,file=paste0(paths$Metadata,"/simulatedTargetParameters.csv"),row.names=FALSE,quote = FALSE,sep=",",col.names=parNames)                #WRITE OPTIMISED PARS TO CSV
#      }
#      write.table(targetMat,file=paste0(paths$Metadata,"/requestedTargets.csv"),row.names=FALSE,quote = FALSE,sep=",")                      #WRITE REQUESTED TARGETS TO CSV
#      write.table(fnam,file = paste0(paths$Scenario,"/filenames.txt"),row.names=FALSE,quote = FALSE,col.names=c("Filename"))                #OUTPUT FILENAME LIST
#    }
#   #----------------------------------------------------------
#
#   #PLOT OF TARGETS TO SUMMARISE
#   if(IOmode == "verbose"){  #If verbose mode
#     #INSERT GENERAL FULL EXPSOURE SPACE PDF
#     # if(!isTRUE(dir.exists(pathPlots))){
#     #   dir.create(pathPlots)
#     # }
#     # if(modelTag[1] != "Simple-ann"){
#     #   exposureSummary(targetMat=simTarget,attSel=attSel,paths=paths)
#     # }else{
#     #   exposureSummary(targetMat=targetMat,attSel=attSel,paths=paths)
#     # }
#   }
#
#   #DIAGNOSTICS GENERATION
#   if(modelTag[1]=="Simple-ann"){ progress("No diagnostics for simple scaling yet",file)}
#
#   #MAKE OUTPUT DATAFRAME IF NOT DEV MODE
#   if(IOmode !="dev"){
#     simDat=list()
#     for(i in 1:nTarget){ # loop over targets
#       simDat[[i]]=makeOutputDataframe(data=sim[[i]],dates=dateExtnd,simVar=simVar,modelTag=modelTag[1])
#     }
#   }
#
#   progress("Creating exposure space diagnostics",file)
#   progress(paste0("Simulations completed. Please see outputs located in ",getwd()),file)
#
#
#    #OUTPUT SIMULATION & ASSOCIATED INFORMATION
#   if(IOmode != "dev"){
#       if(modelTag[[1]] == "Simple-ann"){
#         out=list(sim=sim,
#                  modelTag=modelTag,
#                  target=targetMat,
#                  attPerturb=attPerturb,
#                  attHold=attHold,
#                  attSel=attSel,
#                  attPrim=attPrim,
#                  data=simDat,
#                  nTarget=nTarget,
#                  exSpArgs=exSpArgs)
#       } else {
#         out=list(sim=sim,
#                  modelTag=modelTag,
#                  target=targetMat,
#                  attPerturb=attPerturb,
#                  attHold=attHold,
#                  attSel=attSel,
#                  attPrim=attPrim,
#                  attObs=attObs,
#                  data=simDat,
#                  nTarget=nTarget,
#                  exSpArgs=exSpArgs)
#       }
#
#   }else{
#     out=NULL  #no output in dev mode
#   }
#
#   return(out)
# }



#' Runs a system model and outputs the system performance
#'
#' \code{runSystemModel} uses time series of hydroclimatic variables generated using the function \code{generateScenarios} as input to a systemModel and
#' collates the system performance for all the targets and replicates in the scenarios.
#' @param sim list; a simulation containing the scenarios generated using the function \code{generateScenarios}.
#' @param systemModel a function; The function runs the system model using climate data in a data.frame as input.
#' The function is expected to be created by the user for specific system models. \code{tankWrapper} is an example system model function available in this package.
#' \code{runSystemModel} calls the function \code{systemModel} with two arguments:
#' \itemize{
#' \item {\code{data}} {: data.frame; the climate data in a data frame with columns named \emph{year} \emph{month} \emph{day} \emph{*variable_name1*} \emph{*variable_name2*}. }
#' \item {\code{systemArgs}} {: list; containing the other arguments required by the system model.\code{systemModel} unpack the arguments from the list and uses them as required. }
#' \item {\code{metrics}} {: string vector; containing the names of the performance metrics that the system model returns. It is recommended that the
#' names also contain the units of the metric. See \code{viewTankMetrics()} for examples.}
#' }
#' @param systemArgs a list; containing the input arguments to \code{systemModel}.
#' @param metrics a string vector; the names of the performance metrics the \code{systemModel} function returns.
#' @details The \code{runSystemModel} function code is structured to be simple and may be used as an example to create scripts that use scenarios
#' generated using \code{generateScenarios} to run system models in other programming languages. Type \code{runSystemModel} to view the function code.
#' The function \code{tankWrapper} in this package may be used as an example to create user defined functions for the \code{systemModel} argument.
#' Refer to \code{tankWrapper} to understand how the \code{systemModel} is expected to use \code{systemArgs} and return the calculated performance metrics.
#' The \code{systemModel} function is expected to return a named list of performance metrics. The elements of the vector should correspond to \code{metrics}.
#' @return The function returns a list containing the performance metrics calculated by the \code{systemModel}. Each element of the list corresponds to a performance metric and is named using the \code{metrics} argument.
#' Each element contains performance values calculated at all the target points in the exposure space in a matrix with nrow corresponding to the targets and ncol corresponding to the replicates.
#' @seealso \code{tankWrapper}, \code{generateScenarios}
#' @examples
#' # Example using tankWrapper as the systemModel
#' #=====================================================
#' \dontrun{
#' # create an exposure space
#' attPerturb <- c("P_ann_tot_m", "P_ann_nWet_m")
#' attHold <- c("P_Feb_tot_m", "P_SON_dyWet_m", "P_JJA_avgWSD_m", "P_MAM_tot_m",
#' "P_DJF_avgDSD_m", "Temp_ann_rng_m", "Temp_ann_avg_m")
#' attPerturbType = "regGrid"
#' attPerturbSamp = c(2, 2)
#' attPerturbMin = c(0.9, 0.9)
#' attPerturbMax = c(1.1, 1.1)
#' expSpace <- createExpSpace(attPerturb = attPerturb,
#'                            attPerturbSamp = attPerturbSamp,
#'                            attPerturbMin = attPerturbMin,
#'                            attPerturbMax = attPerturbMax,
#'                            attPerturbType = attPerturbType,
#'                            attHold = attHold,
#'                            attTargetsFile = NULL)
#' # load example observed data available in foreSIGHT
#' data(tankDat)
#' # perform stochastic simulation
#' sim <- generateScenarios(reference = tank_obs,
#'                                   expSpace = expSpace,
#'                                   simLengthNyrs = 30)
#' # use the simulation to run a system model
#' systemArgs <- list(roofArea = 205, nPeople = 1, tankVol = 2400,
#' firstFlush = 2.0, write.file = FALSE)
#' tankMetrics <- viewTankMetrics()
#' systemPerf = runSystemModel(sim = sim,
#'                             systemModel = tankWrapper,
#'                             systemArgs = systemArgs,
#'                             metrics = tankMetrics[1:2])
#'                             }
#' @export


# FUNCTION TO JUST MANAGE THE SIMULATION OF PERFORMANCE
# function may also be named getSystemPerformance
# may use elipsis for unlimited args to systemModel
runSystemModel <- function(sim,                  # output from scenario generator
                           systemModel,          # system model function with arguments
                           systemArgs,           # arguments of the system model
                           metrics               # names of performance metrics returned
){

  # unpacking sim
  repNames <- names(sim[grep("Rep", names(sim))])
  tarNames <- names(sim[[repNames[1]]])

  nRep <- length(repNames)
  nTar <- length(tarNames)

  # variable names (currently fSVars = c("P", "Temp", "PET", "Radn"))
  temp <- names(sim[[repNames[1]]][[tarNames[1]]])
  varNames <- temp[which(temp %in% fSVars)]
  rm(temp)

  # assuming that the system model takes in climate data in a data.frame of the form
  # c(year, month, day, var1, var2)

  performance <- vector("list", length = length(metrics))
  for (i in 1:length(metrics)) performance[[i]] <- matrix(NA, nrow = nTar, ncol = nRep)

  for(r in 1:nRep) {
    for (t in 1:nTar) {
      # initialising scenarioData
      scenarioData <- sim[["simDates"]]
      varTemp <- list()
      for (v in varNames) {
        if ((is.character(sim[["controlFile"]]))) {
          if (sim[["controlFile"]] == "scaling") {
            # using data.frame for single site & list for multi-site
            max_nSites <- max(sapply(sim[[repNames[r]]][[tarNames[t]]], ncol))
            if (max_nSites==1) {
              varTemp <- as.data.frame(sim[[repNames[r]]][[tarNames[t]]][[v]])
            } else{
              varTemp[[v]] <- sim[[repNames[r]]][[tarNames[t]]][[v]]
            }
          } else {
            stop(paste0("sim$controlFile unrecognized."))
          }
        } else {
          # using data.frame for single site & list for multi-site
#          max_nSites <- max(sapply(sim[[repNames[r]]][[tarNames[t]]][[v]], function(x){ncol(x[["sim"]])}))
#          if (max_nSites==1) {
            varTemp <- as.data.frame(sim[[repNames[r]]][[tarNames[t]]][[v]][["sim"]])
#          } else {
#            varTemp[[v]] <- sim[[repNames[r]]][[tarNames[t]]][[v]][["sim"]]
#          }
        }
        if (is.data.frame(varTemp)) {
          names(varTemp) <- v
          scenarioData <- cbind(scenarioData, varTemp)
          rm(varTemp)
        } else{
          scenarioData <- varTemp
        }

      }
      # run the systemModel
      perfTemp <- systemModel(data = scenarioData, systemArgs = systemArgs, metrics = metrics)

      # store performance metrics
      for (i in 1:length(metrics)) performance[[i]][t, r] <- perfTemp[[i]]
      rm(scenarioData)
    }
  }
  names(performance) <- metrics

  return(performance)
}


# function to delete the timeseries in the simulation and keep only the metadata
# this function is intentionally kept different from runSystemModel for simplicity
# may also be named getSimSummary or summarizeSim or getSimMeta or getSimMetaData
# ideally two arguments - 1) a simulation 2) an existing summary metadata object (if NULL existing data will not be added to object)
# if a desired simulation is too large to hold in memory, the function can be used to progressively save metadata into a single summary file

#' Produces a summary object containing the metadata of a full simulation
#'
#' \code{getSimSummary} uses a full simulation generated using the function \code{generateScenarios} as input and outputs a summary object containing the
#' metadata of the full simulation. The output summary object may be used as an input to the plotting functions in this package. The output summary object
#' will be much smaller in size than the full simulation for ease of storage and use with the plotting functions.
#' @param sim list; a simulation containing the scenarios generated using the function \code{generateScenarios}.
#' @seealso \code{generateScenarios}, \code{plotPerformanceSpace}, \code{plotPerformanceOAT}
#' @export

getSimSummary <- function(sim) {

  # if using multiple simulations
  # check that all sim controlFile and expSpace are the same


  # unpacking sim metadata
  repNames <- names(sim[grep("Rep", names(sim))])
  if (is.null(repNames)) stop("There are no replicates in sim.")
  tarNames <- names(sim[[repNames[1]]])
  nRep <- length(repNames)
  nTar <- length(tarNames)

  # variable names
  if (is.list(sim[["controlFile"]])) {
    varNames <- names(sim[["controlFile"]][["modelType"]])
  } else if (!(sim[["controlFile"]]  == "scaling")) {
    stop("sim$controlFile is not recognized.")
  }

  # subsetting sim to simSummary
  simSummary <- list()
  for (r in 1:nRep) {
    if (!is.list(sim[["controlFile"]])) {
      simSummary[[repNames[r]]] <- NULL
    } else {
      for (t in 1:nTar) {
        # save all info except the simulation variables
        indVars <- which(names(sim[[repNames[r]]][[tarNames[t]]]) %in% varNames)
        #othrVars <- names(sim[[repNames[r]]][[tarNames[t]]])[-indVars]
        simSummary[[repNames[r]]][[tarNames[t]]] <- sim[[repNames[r]]][[tarNames[t]]][-indVars]

        # for (i in 1:length(othrVars)) {
        #   simSummary[[repNames[r]]][[tarNames[t]]][[othrVars[i]]] <- sim[[repNames[r]]][[tarNames[t]]][[othrVars[i]]]
        # }
        #
        simSummary[[repNames[r]]][[tarNames[t]]][["seed"]] <- sim[[repNames[r]]][[tarNames[t]]][[varNames[1]]][["seed"]]
      }
    }
  }
  # other metadata
  simSummary[["simDates"]] <- sim[["simDates"]]
  simSummary[["expSpace"]] <- sim[["expSpace"]]
  simSummary[["controlFile"]] <- sim[["controlFile"]]

  return(simSummary)
}






