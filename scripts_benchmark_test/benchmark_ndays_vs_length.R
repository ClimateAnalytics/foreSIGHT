# Passing ndays as an argument versus calculating it as the length of the random vector
#------------------------------------------------------------------------------------------

rm(list = ls())

# Input Arguments
#----------------------
seed <- 1234

ndays <- 3650
sigma <- 5
mu <- -3
lambda <- 1.2
alpha <- 0.8

set.seed(seed)
randomVector <- stats::runif(ndays)

parSigma <- rep_len(sigma, length.out = ndays)
parMu <- rep_len(mu, length.out = ndays)
parLambda <- rep_len(lambda, length.out = ndays)
parAlpha <- rep_len(alpha, length.out = ndays)
parTS <- list(sigma = parSigma, mu = parMu, lambda = parLambda, alpha = parAlpha)


# Function 1
#----------------------
P_latent <- function(parTS,                  # list of parameters
                     randomVector = NULL
) {
  
  # Unpack WGEN parameters
  parAlpha <- parTS$alpha
  parSigma <- parTS$sigma
  parMu <- parTS$mu
  parLambda <- parTS$lambda
  ndays <- length(randomVector)
  
  # Calculate latent variable - latentX
  epsilonT <- stats::qnorm(randomVector, mean = 0, sd = parSigma)
  latentX <- latentX_calc_cpp(parAlpha, epsilonT, ndays)
  latentX <- latentX + parMu
  
  # Transform latentX to rainfall
  rain <- rep_len(0, ndays)  
  latentX_pos_ind <- which(latentX > 0)
  rain[latentX_pos_ind] <- latentX[latentX_pos_ind] ^ parLambda[latentX_pos_ind]
  
  syntP <- list(sim = rain)       
  return(syntP)
  
}

# Function 2
#------------------
P_latent_wndays <- function(parTS,                  # list of parameters
                     randomVector = NULL,
                     ndays
) {
  
  # Unpack WGEN parameters
  parAlpha <- parTS$alpha
  parSigma <- parTS$sigma
  parMu <- parTS$mu
  parLambda <- parTS$lambda
  
  # Calculate latent variable - latentX
  epsilonT <- stats::qnorm(randomVector, mean = 0, sd = parSigma)
  latentX <- latentX_calc_cpp(parAlpha, epsilonT, ndays)
  latentX <- latentX + parMu
  
  # Transform latentX to rainfall
  rain <- rep_len(0, ndays)  
  latentX_pos_ind <- which(latentX > 0)
  rain[latentX_pos_ind] <- latentX[latentX_pos_ind] ^ parLambda[latentX_pos_ind]
  
  syntP <- list(sim = rain)       
  return(syntP)
  
}


ndays_bm <- bench::mark(P_latent(parTS, randomVector), P_latent_wndays(parTS, randomVector, ndays))
ndays_bm

# BM Result: No change in run time
