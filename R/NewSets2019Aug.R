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


GenerateAllRandom200 <- function(parm,
                                 top.level.dir,
                                 mut.mean,
                                 mut.sd,
                                 overwrite = TRUE,
                                 verbose = TRUE) {
  MustCreateDir(top.level.dir, overwrite)
  for (i in 1:nrow(parm)) {
    num.sigs <- parm[i, "total.num.sigs"]
    GenerateOneRowRandom(
      row         = unlist(parm[i, ]),
      dir         = file.path(top.level.dir, paste0(".", num.sigs, ".sigs")),
      num.spectra = 200,
      mut.mean    = mut.mean,
      mut.sd      = mut.sd,
      overwrite   = overwrite,
      verbose     = verbose)
  }
}

GenerateAllRandomSP <-
  function(top.level.dir = "foo", overwite = TRUE, verbose = TRUE) {
  GenerateAllRandom200(
    parm = ParametersSPLike(),
    top.level.dir = top.level.dir,
    mut.mean      = 2.97,
    mut.sd        = 1.331
  )
}


GenerateOneRowRandom <- function(row,
                                 dir,
                                 num.spectra,
                                 mut.mean,
                                 mut.sd,
                                 overwrite,
                                 verbose) {
  MustCreateDir(dir, overwrite)
  total.num.sigs <- row["total.num.sigs"]

  retval <- CreateOneSetOfRandomCatalogs96(
    num.syn.tumors = num.spectra,
             total.num.sigs         = total.num.sigs,
             mut.mean                = mut.mean,
             mut.sd                  = mut.sd,
             mean.num.sigs.per.tumor = row["mean.num.sigs.per.tumor"],
             sd.num.sigs.per.tumor   = row["sd.num.sigs.per.tumor"],
             sig.name.prefix         = "RandSig",
             sample.name.prefix      = "S",
             dir.name                = dir,
             overwrite               = overwrite,
             verbose                 = verbose)

  return(retval)
}


CreateRandomExposures <- function(num.exposures,
                                  mean.num.sigs.per.tumor,
                                  sd.num.sigs.per.tumor,
                                  total.num.sigs,
                                  sig.info,
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
    message("\nCreateOneSetOfRandomCatalogs\n",
            "statistics on number of exposures per tumor")
    message("for targets\nmean(num sigs per tumor) = ",
            mean.num.sigs.per.tumor , "\n",
            "\nsd(num sigs per tumor) = ",
            sd.num.sigs.per.tumor,
            "\ntotal.num.sigs = ", total.num.sigs)
    message("\nNumber of tumors is ", num.exposures)
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
               x, colnames(sigs), sig.info$syn.mean, sig.info$syn.sd) })

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
#' @param mean.num.sigs.per.tumor Mean number of signatures contributing to each tumor.
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

CreateOneSetOfRandomCatalogs96 <-
  function(num.syn.tumors,
           total.num.sigs,
           mut.mean,
           mut.sd,
           mean.num.sigs.per.tumor,
           sd.num.sigs.per.tumor,
           sig.name.prefix,
           sample.name.prefix,
           dir.name,
           overwrite = FALSE,
           verbose = TRUE) {

    sigs <-
      CreateRandomMutSigProfiles(
        ICAMS::catalog.row.order[["SBS96"]], total.num.sigs, sig.name.prefix)

    sig.info <- CreateMeanAndStdevForSigs(
      total.num.sigs, mut.mean, mut.sd, colnames(sigs))

    exp <-
      CreateRandomExposures(
        num.exposures           = num.syn.tumors,
        mean.num.sigs.per.tumor = mean.num.sigs.per.tumor,
        sd.num.sigs.per.tumor   = sd.num.sigs.per.tumor,
        total.num.sigs          = total.num.sigs,
        sig.info                = sig.info,
        sample.name.prefix      = sample.name.prefix,
        sigs                    = sigs,
        verbose                 = verbose)

    NewCreateAndWriteCatalog(
      sigs      = sigs,
      exp       = exp,
      dir       = dir.name,
      overwrite = overwrite)

  }
