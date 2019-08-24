#' Parameters for generating data with exposures resembling SigProfile exposures.
#' @keywords internal

ParametersSPLike <- function() {
  retval <- data.frame(
    total.num.sigs          = c(3, 5, 10, 15, 25),

    # Except for the first element, these are based on the observed mean and
    # standard deviation in PCAWG tumors.
    mean.num.sigs.per.tumor = c(1.5, rep(3.9, 4)),
    sd.num.sigs.per.tumor   = c(0.5, rep(1.3, 4))
  )
  return(retval)
}


#' Parameters for generating data with exposures resembling SignatureAnalyzer exposures.
#' @keywords internal

ParametersSALike <- function() {
  retval <- data.frame(
    total.num.sigs          = c(3, 5, 10, 15, 25),

    # These are chosen to ramp up to the mean standard deviation in
    # PCAWG tumors (15.525 and 6.172).
    mean.num.sigs.per.tumor = c(2.9, 4.9, 8, 11, 16),
    sd.num.sigs.per.tumor   = c(1,   2,   4,  5,  6))
  return(retval)
}


CreateRandomExposures <- function(num.exposures,
                                  mean.num.sigs.per.tumor,
                                  sd.num.sigs.per.tumor,
                                  total.num.sigs,
                                  per.sig.mean.and.sd,
                                  sample.name.prefix,
                                  sigs,
                                  verbose) {

  buffer <- 100

  exp.nums <-
    CreateExposuresNums(
      num.exposures = num.exposures + buffer,
      mean = mean.num.sigs.per.tumor,
      sd = sd.num.sigs.per.tumor,
      total.num.sigs = total.num.sigs)

  if (verbose) {
      message("\n", paste(rep("=", 40), collapse=""),
              "\nCreate1CatRandomExpKnownSigs",
              "\ntotal.num.sigs          = ", total.num.sigs,
              "\nmean.num.sigs.per.tumor = ", mean.num.sigs.per.tumor,
              "\nsd.num.sigs.per.tumor   = ", sd.num.sigs.per.tumor)
    ss <- summary(exp.nums)
    for (nn in names(ss)) {
      message(nn, " = ", ss[nn])
    }
    message("sd = ", sd(exp.nums), "\n")
  }

  exp <-
    sapply(exp.nums,
           function(x) {
             ExposureNums2Exposures(
               x,
               colnames(sigs),
               per.sig.mean.and.sd$syn.mean,
               per.sig.mean.and.sd$syn.sd)
           })

  test.catalog <- sigs %*% exp
  test.catalog <- round(test.catalog, digits = 0)
  zero.mutations <- colSums(test.catalog) == 0
  # colSums(test.catalog) == 0 can occur after rounding even if
  # any(colSUms(test.catalog) < 1) before rounding is FALSE, if
  # before rounding mutiple mutational classes had < 0.5 mutations.


  exp <- exp[ , !zero.mutations]
  if (ncol(exp) < num.exposures)
    stop("Too many tumors with no mutations; check the code, ",
         "possibly increase the value of variable buffer")
  exp <- exp[ , 1:num.exposures]
  colnames(exp) <- paste0(sample.name.prefix, 1:num.exposures)

  # stopifnot(!any(colSums(test.catalog) < 1))

  return(exp)
}



GenerateAllRandom200 <- function(parm,
                                 top.level.dir,
                                 mut.mean,
                                 mut.sd,
                                 overwrite = TRUE,
                                 verbose = TRUE) {
  MustCreateDir(top.level.dir, overwrite)
  for (i in 1:nrow(parm)) {
    num.sigs <- parm[i, "total.num.sigs"]
    # Replace with Create1CatRandomExpRandomSigs Todo
    Generate1RowRandom(
      row         = unlist(parm[i, ]),
      dir.name    = file.path(top.level.dir, paste0(".", num.sigs, ".sigs")),
      num.spectra = 200,
      mut.mean    = mut.mean,
      mut.sd      = mut.sd,
      overwrite   = overwrite,
      verbose     = verbose)
  }
}

GenerateAllRandomSA <-
  function(top.level.dir = "foo", overwite = TRUE, verbose = TRUE) {
    GenerateAllRandom200(
      parm = ParametersSALike(),
      top.level.dir = top.level.dir,
      mut.mean      = 2.349, # Change mean.log10.mut.per.sig
      # Based on mean(log10(sa.all.real.exposures[sa.all.real.exposures >= 1]))

      mut.sd        =  0.6641
      # Based on sd(log10(sa.all.real.exposures[sa.all.real.exposures >= 1]))
      )
  }



GenerateAllRandomSP <-
  function(top.level.dir = "foo", overwite = TRUE, verbose = TRUE) {
    GenerateAllRandom200(
      parm = ParametersSPLike(),
      top.level.dir = top.level.dir,
      mut.mean      = 2.97,
      # Based on mean(log10(sp.all.real.exposures[sp.all.real.exposures >= 1]))

      mut.sd        = 0.7047
      # Based on sd(log10(sp.all.real.exposures[sp.all.real.exposures >= 1]))
    )
  }


Generate1RowRandom <- function(row,
                               num.spectra,
                               mut.mean,
                               mut.sd,
                               overwrite,
                               verbose,
                               dir.name) {

  retval <- Create1CatRandomExpRandomSigs(
    num.syn.tumors = num.spectra,
    total.num.sigs         = row["total.num.sigs"],
    mut.mean                = mut.mean,
    mut.sd                  = mut.sd,
    mean.num.sigs.per.tumor = row["mean.num.sigs.per.tumor"],
    sd.num.sigs.per.tumor   = row["sd.num.sigs.per.tumor"],
    dir.name                = dir.name,
    overwrite               = overwrite,
    verbose                 = verbose)
  return(retval)
}


#' Create catalog of "random" synthetic spectra for 96-channel mutation types.
#'
#' @param num.syn.tumors Total number of synthetic tumors to create.
#'
#' @param total.num.sigs Total number of signatures in the universe.
#'
#' @param mut.mean Mean of the log10 of the
#' number of mutations due to each signature.
#'
#' @param mut.sd Standard deviation of the log10 of
#'  the number of mutations due to each signature.
#'
#' @param mean.num.sigs.per.tumor Mean number of signatures
#'        contributing to each tumor.
#'
#' @param sd.num.sigs.per.tumor Standard deviation the number of signatures
#' in each tumor.
#'
#' @param sig.name.prefix String to put in front of an integer (as
#' string) to form an identifier for a synthetic signature.
#'
#' @param sample.name.prefix String to put in front of an integer (as
#' string) to form an identifier for a synthetic sample (tumor).
#'
#' @param dir.name A string indicating the name of the 96-channel
#' subdirectory; probably one of \code{"sa.sa.COMPOSITE"} or
#' \code{"sp.sa.COMPOSITE"}.
#'
#' @param COMPOSITE.features Character vector containing
#' rownames for a COMPOSITE signature or catalog.
#'
#' @param overwrite If \code{TRUE} overwrite existing directories / files.
#'
#' @keywords internal

Create1CatRandomExpRandomSigs <-
  function(num.syn.tumors,
           total.num.sigs,
           mut.mean,
           mut.sd,
           mean.num.sigs.per.tumor,
           sd.num.sigs.per.tumor,
           dir.name,
           overwrite = FALSE,
           verbose = TRUE) {

    sigs <-
      CreateRandomMutSigProfiles(
        ICAMS::catalog.row.order[["SBS96"]], total.num.sigs,
        sig.name.prefix = "RandSig")

    sig.info <- CreateMeanAndStdevForSigs(
      total.num.sigs, mut.mean, mut.sd, colnames(sigs))

    exp <-
      CreateRandomExposures(
        num.exposures           = num.syn.tumors,
        mean.num.sigs.per.tumor = mean.num.sigs.per.tumor,
        sd.num.sigs.per.tumor   = sd.num.sigs.per.tumor,
        total.num.sigs          = total.num.sigs,
        per.sig.mean.and.sd     = sig.info,
        sample.name.prefix      = "S",
        sigs                    = sigs,
        verbose                 = verbose)

    NewCreateAndWriteCatalog(
      sigs      = sigs,
      exp       = exp,
      dir       = dir.name,
      overwrite = overwrite)

  }

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
  logfile <- file.path(top.level.dir, "log.txt")

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

