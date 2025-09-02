# unit tests limited to annual LV model.

test_that("SWGsim.latent output matches saved reference", {
  randomVector <- c(0.5349710, 0.5463762, 0.8611103, 0.6066115, 0.6701731)
  randomUnitNormalVector <- stats::qnorm(randomVector)
  parSigma <- rep_len(5, length.out = 5)
  parMu <- rep_len(-3, length.out = 5)
  parLambda <- rep_len(1.2, length.out = 5)
  parAlpha <- rep_len(0.8, length.out = 5)

  parTS <- list(sigma = parSigma, mu = parMu, lambda = parLambda, alpha = parAlpha)

  actual <- SWGsim.latent(SWGpar = parTS, randomTerm = list(randomUnitNormalVector = randomUnitNormalVector), nTimes = 5)

  expected <- readRDS(testthat::test_path("../P_latent_output1.rds"))$sim

  expect_equal(actual, expected)
})

test_that("simClim output matches saved reference", {
  modelTag <- "P-ann-latent"
  modelInfo <- modelInfoList[[modelTag]]

  datInd <- list()
  datInd$nTimes <- 5

  parS <- c(0.8, 5, -3, 1.2)
  randomVector <- c(0.5349710, 0.5463762, 0.8611103, 0.6066115, 0.6701731)

  actual <- simClim(
    parS = parS,
    modelTag = modelTag,
    modelInfo = modelInfo,
    datInd = datInd,
    randomTerm = list(randomVector = randomVector)
  )

  expected <- readRDS(testthat::test_path("../P_latent_output1.rds"))$sim

  expect_equal(actual, expected)
})
