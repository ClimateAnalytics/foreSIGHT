# Inverting random vector (q****) versus random number generation (r****)
#--------------------------------------------------------------------------

# 1. Normal
#------------------

rm(list = ls())

seed <- 5312
parSigma <- 7.1
ndays <- 9000

f1 <- function() {
  set.seed(seed)
  stats::rnorm(n = ndays, mean = 0, sd = parSigma)
}

randomVector <- pnorm(f1(), mean = 0, sd = parSigma)
f2 <- function() { stats::qnorm(randomVector, mean = 0, sd = parSigma) }

norm_bm <- bench::mark(f1(), f2())
norm_bm

# BM Results
# expression      min   median
# f1()        184.2µs    204µs
# f2()         79.3µs   84.6µs

# bm result: Inversion of normal distribution is about 2.2 times faster than normal random number generation


# 2. Gamma
#-------------------

rm(list = ls())

seed <- 7613
parShape <- 2.3
parScale <- 4.1
ndays <- 3650

f1 <- function() {
  set.seed(seed)
  rgamma(n = ndays, shape = parShape, scale = parScale)
}

randomVector <- pgamma(f1(), shape = parShape, scale = parScale)
f2 <- function() { qgamma(randomVector, shape = parShape, scale = parScale) }

gamma_bm <- bench::mark(f1(), f2())
gamma_bm

# BM Results
# expression      min   median
# f1()       308.85µs 320.12µs
# f2()         4.29ms   4.77ms 

# bm result: Gamma inversion is about 14-15 times slower than gamma random number generation








