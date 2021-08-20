context("Generate synthetic tumors from a list of signature parameters")

test_that("Generate synthetic exposures for SBS96", {
  input.sigs.SBS96 <- PCAWG7::signature$genome$SBS96
  real.exposures.SBS96 <- PCAWG7::exposure$PCAWG$SBS96
  cancer.types <- PCAWG7::CancerTypes()[1:5]
  sig.params <- SynSigGen::signature.params$SBS96
  SBS96.sig.params <-
    GenerateListOfSigParams(real.exposures = real.exposures.SBS96,
                            cancer.types = cancer.types,
                            distribution = "neg.binom",
                            sig.params = sig.params
    )
  dirs <- file.path(tempdir(), paste0("SBS96.synthetic.tumors", 1:2))

  SBS96.synthetic.tumors1 <-
    GenerateSyntheticTumorsFromSigParams(seed = 191906,
                                         dir = dirs[1],
                                         cancer.types = cancer.types,
                                         samples.per.cancer.type = 30,
                                         input.sigs = input.sigs.SBS96,
                                         sig.params = SBS96.sig.params,
                                         distribution = "neg.binom",
                                         sample.prefix.name = "SP.Syn."
    )

  SBS96.synthetic.tumors2 <-
    GenerateSyntheticTumors(seed = 191906,
                            dir = dirs[2],
                            cancer.types = cancer.types,
                            samples.per.cancer.type = 30,
                            input.sigs = input.sigs.SBS96,
                            real.exposures = real.exposures.SBS96,
                            distribution = "neg.binom",
                            sample.prefix.name = "SP.Syn."
    )

  expect_equal(SBS96.synthetic.tumors1$ground.truth.exposures,
               SBS96.synthetic.tumors2$ground.truth.exposures)

  sapply(dirs, FUN = unlink, recursive = TRUE)
}
)

test_that("Generate synthetic exposures for DBS78", {
  input.sigs.DBS78 <- PCAWG7::signature$genome$DBS78
  real.exposures.DBS78 <- PCAWG7::exposure$PCAWG$DBS78
  cancer.types <- PCAWG7::CancerTypes()[1:5]
  sig.params <- SynSigGen::signature.params$DBS78
  DBS78.sig.params <-
    GenerateListOfSigParams(real.exposures = real.exposures.DBS78,
                            cancer.types = cancer.types,
                            distribution = "neg.binom",
                            sig.params = sig.params
    )
  dirs <- file.path(tempdir(), paste0("DBS78.synthetic.tumors", 1:2))

  DBS78.synthetic.tumors1 <-
    GenerateSyntheticTumorsFromSigParams(seed = 191906,
                                         dir = dirs[1],
                                         cancer.types = cancer.types,
                                         samples.per.cancer.type = 30,
                                         input.sigs = input.sigs.DBS78,
                                         sig.params = DBS78.sig.params,
                                         distribution = "neg.binom",
                                         sample.prefix.name = "SP.Syn."
    )

  DBS78.synthetic.tumors2 <-
    GenerateSyntheticTumors(seed = 191906,
                            dir = dirs[2],
                            cancer.types = cancer.types,
                            samples.per.cancer.type = 30,
                            input.sigs = input.sigs.DBS78,
                            real.exposures = real.exposures.DBS78,
                            distribution = "neg.binom",
                            sample.prefix.name = "SP.Syn."
    )

  expect_equal(DBS78.synthetic.tumors1$ground.truth.exposures,
               DBS78.synthetic.tumors2$ground.truth.exposures)

  sapply(dirs, FUN = unlink, recursive = TRUE)
}
)

test_that("Generate synthetic exposures for ID", {
  input.sigs.ID <- PCAWG7::signature$genome$ID
  real.exposures.ID <- PCAWG7::exposure$PCAWG$ID
  cancer.types <- PCAWG7::CancerTypes()[1:5]
  sig.params <- SynSigGen::signature.params$ID
  ID.sig.params <-
    GenerateListOfSigParams(real.exposures = real.exposures.ID,
                            cancer.types = cancer.types,
                            distribution = "neg.binom",
                            sig.params = sig.params
    )
  dirs <- file.path(tempdir(), paste0("ID.synthetic.tumors", 1:2))

  ID.synthetic.tumors1 <-
    GenerateSyntheticTumorsFromSigParams(seed = 191906,
                                         dir = dirs[1],
                                         cancer.types = cancer.types,
                                         samples.per.cancer.type = 30,
                                         input.sigs = input.sigs.ID,
                                         sig.params = ID.sig.params,
                                         distribution = "neg.binom",
                                         sample.prefix.name = "SP.Syn."
    )

  ID.synthetic.tumors2 <-
    GenerateSyntheticTumors(seed = 191906,
                            dir = dirs[2],
                            cancer.types = cancer.types,
                            samples.per.cancer.type = 30,
                            input.sigs = input.sigs.ID,
                            real.exposures = real.exposures.ID,
                            distribution = "neg.binom",
                            sample.prefix.name = "SP.Syn."
    )

  expect_equal(ID.synthetic.tumors1$ground.truth.exposures,
               ID.synthetic.tumors2$ground.truth.exposures)

  sapply(dirs, FUN = unlink, recursive = TRUE)
}
)
