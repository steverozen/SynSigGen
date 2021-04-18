context("Generate synthetic exposures for DBS78 signatures")

test_that("Generate synthetic exposures for DBS78 signatures", {
  real.exposures.DBS78 <- PCAWG7::exposure$PCAWG$DBS78
  params.DBS78 <- GetSynSigParamsFromExposures(exposures = real.exposures.DBS78)
  synthetic.exposures.DBS78 <- GenerateSyntheticExposures(sig.params = params.DBS78)
  sigs.DBS78 <- PCAWG7::signature$genome$DBS78
  out <- CreateAndWriteCatalog(sigs = sigs.DBS78,
                               exp = synthetic.exposures.DBS78,
                               my.dir = file.path(tempdir(), "DBS78"))

}
)
