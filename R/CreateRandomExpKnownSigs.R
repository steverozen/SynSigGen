# ============================================================
# Generate spectra from randomly selected sets of
# real mutational signatures.

MeanAndSD1Sig <- function(exposure.row) {
  exp <- exposure.row[exposure.row > 0.5]
  exp <- log10(exp)
  return(c(mean = mean(exp), sd = sd(exp)))
}


PerSigMeanAndSDSA <- function()  {
  res1 <- apply(SynSigGen::sa.all.real.exposures, MARGIN = 1, FUN = MeanAndSD1Sig)
  res1 <- res1[ , !is.na(res1[2, ])]
  return(res1)
}


PerSigMeanAndSDSP <- function() {
  res1 <- apply(SynSigGen::sp.all.real.exposures, MARGIN = 1, FUN = MeanAndSD1Sig)
  res1 <- res1[ , !is.na(res1[2, ])]
  return(res1)
}


Generate.SA.signatures.random.subsets <-
  function(top.level.dir = "data-raw/SA.signatures.random.subsets",
           overwrite     = TRUE,
           verbose       = TRUE,
           seed          = 3100) {
    MustCreateDir(top.level.dir, overwrite = TRUE)
    log <- testthat::capture_messages(
      GenerateMatrixRandomExpKnownSigs200(
        parm          = ParametersSALike(),
        top.level.dir = top.level.dir,
        sigs          = SynSigGen::sa.96.sigs,
        sig.info      = PerSigMeanAndSDSA(),
        overwrite     = overwrite,
        verbose       = verbose,
        seed          = seed
      ))
    write(x = log, file = file.path(top.level.dir, "log.txt"), append = TRUE)
  }

Generate.SP.signatures.random.subsets <-
  function(top.level.dir = "data-raw/SP.signatures.random.subsets",
           overwrite     = TRUE,
           verbose       = TRUE,
           seed          = 3120) {
    MustCreateDir(top.level.dir, overwrite = TRUE)
    log <- testthat::capture_messages(
      GenerateMatrixRandomExpKnownSigs200(
        parm          = ParametersSPLike(),
        top.level.dir = top.level.dir,
        sigs          = SynSigGen::sp.sigs,
        sig.info      = PerSigMeanAndSDSP(),
        overwrite     = overwrite,
        verbose       = verbose,
        seed          = seed))
    write(x = log, file = file.path(top.level.dir, "log.txt"), append = TRUE)

  }


GenerateMatrixRandomExpKnownSigs200 <- function(parm,
                                                top.level.dir,
                                                overwrite = TRUE,
                                                verbose = TRUE,
                                                sigs,
                                                sig.info,
                                                seed = seed) {
  MustCreateDir(top.level.dir, overwrite)

  # RNGMessages("GenerateMatrixRandomExpKnownSigs200 before")
  set.seed(seed, sample.kind = "Rejection")
  RNGMessages("GenerateMatrixRandomExpKnownSigs200 after", cat)

  for (i in 1:nrow(parm)) {
    Create1CatRandomExpKnownSigs(
      top.level.dir = top.level.dir,
      num.syn.tumors = 200,
      parm.row       = parm[i, ],
      sigs           = sigs,
      sig.info       = sig.info,
      overwrite      = overwrite,
      verbose        = verbose)
  }
}


Create1CatRandomExpKnownSigs <-
  function(top.level.dir,
           num.syn.tumors,
           parm.row,
           sigs,
           sig.info,
           overwrite = FALSE,
           verbose = TRUE) {

    parm.row <- unlist(parm.row)
    total.num.sigs <- parm.row["total.num.sigs"]

    dir <-  file.path(top.level.dir, paste0(total.num.sigs, ".sigs"))

    MustCreateDir(dir, overwrite)

    # Select total.num.sigs columns from sig.info
    cols.to.use <- sample(ncol(sig.info), total.num.sigs)
    sig.info    <- sig.info[ , cols.to.use, drop = FALSE]

    per.sig.mean.and.sd <- list(syn.mean = sig.info[1, ],
                                syn.sd   = sig.info[2, ])

    sigs <- sigs[ , colnames(sig.info)] # We do not have sig info
                                        # for some rare signatures.
                                        # These rare signatures are not
                                        # in sig.info.

    exp <-
      CreateRandomExposures(
        num.exposures           = num.syn.tumors,
        mean.num.sigs.per.tumor = parm.row["mean.num.sigs.per.tumor"],
        sd.num.sigs.per.tumor   = parm.row["sd.num.sigs.per.tumor"],
        total.num.sigs          = total.num.sigs,
        per.sig.mean.and.sd     = per.sig.mean.and.sd,
        sample.name.prefix      = "S",
        sigs                    = sigs,
        verbose                 = verbose)

    NewCreateAndWriteCatalog(
      sigs      = sigs,
      exp       = exp,
      dir       = dir,
      overwrite = overwrite)

  }

