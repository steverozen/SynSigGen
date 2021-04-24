context("Generate synthetic tumors from different distributions")

test_that("Generate synthetic exposures for DBS78 signatures using log-normal distribution", {
  set.seed(2586)
  real.exposures.DBS78 <- PCAWG7::exposure$PCAWG$DBS78
  params.DBS78 <- GetSynSigParamsFromExposures(exposures = real.exposures.DBS78)
  synthetic.exposures.DBS78 <- GenerateSyntheticExposures(sig.params = params.DBS78)
  sigs.DBS78 <- PCAWG7::signature$genome$DBS78
  out <- CreateAndWriteCatalog(sigs = sigs.DBS78,
                               exp = synthetic.exposures.DBS78,
                               my.dir = file.path(tempdir(), "DBS78"),
                               overwrite = TRUE)

}
)

test_that("Generate synthetic exposures for DBS78 signatures using negative binomial distribution", {
  set.seed(2586)
  real.exposures.DBS78 <- PCAWG7::exposure$PCAWG$DBS78
  sigs.DBS78 <- PCAWG7::signature$genome$DBS78
  params.DBS78.neg.binom <-
    GetSynSigParamsFromExposures(exposures = real.exposures.DBS78,
                                 distribution = "neg.binom")
  synthetic.exposures.DBS78.neg.binom <-
    GenerateSyntheticExposures(sig.params = params.DBS78.neg.binom,
                               distribution = "neg.binom")
  out <- CreateAndWriteCatalog(sigs = sigs.DBS78,
                               exp = synthetic.exposures.DBS78.neg.binom,
                               my.dir = file.path(tempdir(), "DBS78.neg.binom"),
                               overwrite = TRUE)

}
)
