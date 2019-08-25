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
        parm           = ParametersSALike(),
        top.level.dir  = top.level.dir,
        sigs           = SynSigGen::sa.96.sigs,
        sig.info       = PerSigMeanAndSDSA(),
        overwrite      = overwrite,
        verbose        = verbose,
        seed           = seed,
        num.replicates = 10
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
        parm           = ParametersSPLike(),
        top.level.dir  = top.level.dir,
        sigs           = SynSigGen::sp.sigs,
        sig.info       = PerSigMeanAndSDSP(),
        overwrite      = overwrite,
        verbose        = verbose,
        seed           = seed,
        num.replicates = 10))
    write(x = log, file = file.path(top.level.dir, "log.txt"), append = TRUE)

  }


GenerateMatrixRandomExpKnownSigs200 <- function(parm,
                                                top.level.dir,
                                                overwrite = TRUE,
                                                verbose = TRUE,
                                                sigs,
                                                sig.info,
                                                seed = seed,
                                                num.replicates = 1) {
  MustCreateDir(top.level.dir, overwrite)

  # RNGMessages("GenerateMatrixRandomExpKnownSigs200 before")
  set.seed(seed, sample.kind = "Rejection")
  RNGMessages("GenerateMatrixRandomExpKnownSigs200 after", cat)

  summary.file <- file.path(top.level.dir, "summary.csv")

  for (i in 1:nrow(parm)) {

    total.num.sigs          <- parm[i, "total.num.sigs"]
    mean.num.sigs.per.tumor <- parm[i, "mean.num.sigs.per.tumor"]
    sd.num.sigs.per.tumor   <- parm[i, "sd.num.sigs.per.tumor"]
    if (verbose) {
      message("\n", paste(rep("=", 40), collapse = ""))
      message("\nGenerateRandomExpRandomSigs200\ntarget distribution parameters:",
              "\ntotal.num.sigs          = ", total.num.sigs,
              "\nmean.num.sigs.per.tumor = ", mean.num.sigs.per.tumor,
              "\nsd.num.sigs.per.tumor   = ", sd.num.sigs.per.tumor)
    }

    for (replicate.number in 1:num.replicates) {

      dir.name <-
        file.path(top.level.dir,
                  paste0("syn.data.", total.num.sigs, ".sigs"))
      if (num.replicates > 1) {
        dir.name <- paste0(dir.name, ".repnum.",
                           formatC(replicate.number, width = 2, flag = "0"))
      }

      if (verbose) {
        message("\n", paste(rep("-", 40), collapse = ""))
        message("replicate number ", replicate.number)
        message("directory is ", dir.name)
      }

      actual.sig.num.mean.and.sd <-
        Create1CatRandomExpKnownSigs(
          num.syn.tumors          = 200,
          total.num.sigs          = total.num.sigs,
          mean.num.sigs.per.tumor = mean.num.sigs.per.tumor,
          sd.num.sigs.per.tumor   = sd.num.sigs.per.tumor,
          dir.name                = dir.name,
          sigs                    = sigs,
          sig.info                = sig.info,
          overwrite               = overwrite,
          verbose                 = verbose)

      cat(total.num.sigs, mean.num.sigs.per.tumor, sd.num.sigs.per.tumor,
          replicate.number, "mumble mean", "mumble sd",
          "\n", sep = ",", file = summary.file, append = TRUE)
    }
  }
}


Create1CatRandomExpKnownSigs <-
  function(num.syn.tumors,
           total.num.sigs,
           mean.num.sigs.per.tumor,
           sd.num.sigs.per.tumor,
           dir.name,
           sigs,
           sig.info,
           overwrite = FALSE,
           verbose = TRUE) {

    parm.row <- unlist(parm.row)

    MustCreateDir(dir.name, overwrite)

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
        mean.num.sigs.per.tumor = mean.num.sigs.per.tumor,
        sd.num.sigs.per.tumor   = sd.num.sigs.per.tumor,
        total.num.sigs          = total.num.sigs,
        per.sig.mean.and.sd     = per.sig.mean.and.sd, # Mean and sd of log10 of distribution of number of mutations due to each signature.
        sample.name.prefix      = "S",
        sigs                    = sigs,
        verbose                 = verbose)

    NewCreateAndWriteCatalog(
      sigs      = sigs,
      exp       = exp,
      dir       = dir.name,
      overwrite = overwrite)

    return(c(
      actual.sig.num.mean = attr(exp, "actual.sig.num.mean"),
      actual.sig.num.sd = attr(exp, "actual.sig.num.sd")))

  }

