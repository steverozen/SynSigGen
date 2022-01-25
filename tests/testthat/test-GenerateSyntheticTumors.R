context("Generate synthetic tumors from different distributions")

test_that("Generate synthetic exposures for DBS78 using log-normal distribution", {
  set.seed(2586)
  real.exposures.DBS78 <- PCAWG7::exposure$PCAWG$DBS78
  params.DBS78 <- GetSynSigParamsFromExposures(exposures = real.exposures.DBS78)
  synthetic.exposures.DBS78 <- GenerateSyntheticExposures(sig.params = params.DBS78)
  sigs.DBS78 <- cosmicsig::COSMIC_v3.2$signature$GRCh37$DBS78
  out <- CreateAndWriteCatalog(sigs = sigs.DBS78,
                               exp = synthetic.exposures.DBS78,
                               my.dir = file.path(tempdir(), "DBS78"),
                               overwrite = TRUE)
  expect_equal(dim(out), c(78, 10))
  unlink(file.path(tempdir(), "DBS78"), recursive = TRUE)

}
)

test_that("Generate synthetic exposures for DBS78 using negative binomial distribution", {
  set.seed(2586)
  real.exposures.DBS78 <- PCAWG7::exposure$PCAWG$DBS78
  sigs.DBS78 <- cosmicsig::COSMIC_v3.2$signature$GRCh37$DBS78
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
  expect_equal(dim(out), c(78, 10))
  unlink(file.path(tempdir(), "DBS78.neg.binom"), recursive = TRUE)

}
)

test_that("Generate synthetic tumors for DBS78 by cancer types using log-normal distribution", {
  # Generate synthetic tumors for DBS78
  input.sigs.DBS78 <- cosmicsig::COSMIC_v3.2$signature$GRCh37$DBS78
  real.exposures.DBS78 <- PCAWG7::exposure$PCAWG$DBS78
  cancer.types <- PCAWG7::CancerTypes()[1:5]
  DBS78.synthetic.tumors <-
    GenerateSyntheticTumors(seed = 191906,
                            dir = file.path(tempdir(), "DBS78.synthetic.tumors"),
                            cancer.types = cancer.types,
                            samples.per.cancer.type = 30,
                            input.sigs = input.sigs.DBS78,
                            real.exposures = real.exposures.DBS78,
                            sample.prefix.name = "SP.Syn."
    )
  expect_equal(dim(DBS78.synthetic.tumors$ground.truth.catalog),
               c(78, 150))
  unlink(file.path(tempdir(), "DBS78.synthetic.tumors"), recursive = TRUE)

}
)

test_that("Generate synthetic tumors for DBS78 by cancer types using negative binomial distribution", {
  # Generate synthetic tumors for DBS78
  input.sigs.DBS78 <- cosmicsig::COSMIC_v3.2$signature$GRCh37$DBS78
  real.exposures.DBS78 <- PCAWG7::exposure$PCAWG$DBS78
  cancer.types <- PCAWG7::CancerTypes()[1:5]
  DBS78.synthetic.tumors <-
    GenerateSyntheticTumors(seed = 191906,
                            dir = file.path(tempdir(), "DBS78.synthetic.tumors.neg.binom"),
                            cancer.types = cancer.types,
                            samples.per.cancer.type = 30,
                            input.sigs = input.sigs.DBS78,
                            real.exposures = real.exposures.DBS78,
                            distribution = "neg.binom",
                            sample.prefix.name = "SP.Syn."
    )
  expect_equal(dim(DBS78.synthetic.tumors$ground.truth.catalog),
               c(78, 150))
  unlink(file.path(tempdir(), "DBS78.synthetic.tumors.neg.binom"), recursive = TRUE)

}
)

test_that("Generate synthetic tumors for ID by cancer types using log-normal distribution", {
  # Generate synthetic tumors for ID
  input.sigs.ID <- cosmicsig::COSMIC_v3.2$signature$GRCh37$ID
  real.exposures.ID <- PCAWG7::exposure$PCAWG$ID
  cancer.types <- PCAWG7::CancerTypes()[1:5]
  ID.synthetic.tumors <-
    GenerateSyntheticTumors(seed = 191906,
                            dir = file.path(tempdir(), "ID.synthetic.tumors"),
                            cancer.types = cancer.types,
                            samples.per.cancer.type = 30,
                            input.sigs = input.sigs.ID,
                            real.exposures = real.exposures.ID,
                            sample.prefix.name = "SP.Syn."
    )
  expect_equal(dim(ID.synthetic.tumors$ground.truth.catalog),
               c(83, 150))
  unlink(file.path(tempdir(), "ID.synthetic.tumors"), recursive = TRUE)

}
)

test_that("Generate synthetic tumors for ID by cancer types using negative binomial distribution", {
  # Generate synthetic tumors for ID
  input.sigs.ID <- cosmicsig::COSMIC_v3.2$signature$GRCh37$ID
  real.exposures.ID <- PCAWG7::exposure$PCAWG$ID
  cancer.types <- PCAWG7::CancerTypes()[1:5]
  ID.synthetic.tumors <-
    GenerateSyntheticTumors(seed = 191906,
                            dir = file.path(tempdir(), "ID.synthetic.tumors.neg.binom"),
                            cancer.types = cancer.types,
                            samples.per.cancer.type = 30,
                            input.sigs = input.sigs.ID,
                            real.exposures = real.exposures.ID,
                            distribution = "neg.binom",
                            sample.prefix.name = "SP.Syn."
    )
  expect_equal(dim(ID.synthetic.tumors$ground.truth.catalog),
               c(83, 150))
  unlink(file.path(tempdir(), "ID.synthetic.tumors.neg.binom"), recursive = TRUE)

}
)
