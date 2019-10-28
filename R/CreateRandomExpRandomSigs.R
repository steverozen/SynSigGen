# CreateRanomExpRandomSigs.R
#
# Generate a suite of data with random signatures and exposures.


#' Parameters for generating data with exposures resembling SignatureAnalyzer exposures.
#' @keywords internal

ParametersSALike <- function() {
  retval <- data.frame(
    total.num.sigs          = c(3,   5,   10,  15,  25),

    # These are chosen to ramp up to the mean standard deviation in
    # PCAWG tumors (15.525 and 6.172).
    mean.num.sigs.per.tumor = c(2.9, 4.9, 8,   11,  16),
    sd.num.sigs.per.tumor   = c(1,   2,   4,    5,   6),
    mean.log10.exp.per.sig  = c(2.6, 2.4, 2.1,  2.1, 2),
    sd.log10.exp.per.sig    = c(0.6, 0.6, 0.5,  0.4, 0.4)
  )


  return(retval)
}


#' Generate a full set of random data with characteristics somewhat like SignatureAnalyzer signatures and attribution.
#' @keywords internal
#'
#' Base parameters on SignatureAnalyzer attributions (exposures).

GenerateAllRandomSA <-
  function(top.level.dir = "../SA.like.rand.exp.rand.sigs.2019.10.28",
             # "../SA.like.rand.exp.rand.sigs.2019.08.26",
           overwite = TRUE,
           verbose = TRUE,
           num.replicates = 10) {
    log <- testthat::capture_messages(
      GenerateRandomExpRandomSigs200(
        parm                   = ParametersSALike(),
        top.level.dir          = top.level.dir,

        # mean.log10.exp.per.sig = 2.349, # Change mean.log10.mut.per.sig
        # Based on mean(log10(sa.all.real.exposures[sa.all.real.exposures >= 1]))

        # sd.log10.exp.per.sig  =  0.6641,
        # Based on sd(log10(sa.all.real.exposures[sa.all.real.exposures >= 1]))

        seed                  = 811211,
        num.replicates        = num.replicates))
    cat(log, file = file.path(top.level.dir, "log.txt"))
    return(Resummarize(top.level.dir))
  }


#' Parameters for generating data with exposures resembling SigProfiler exposures.
#' @keywords internal

ParametersSPLike <- function() {
  retval <- data.frame(
    total.num.sigs          = c(3, 5, 10, 15, 25),

    # Except for the first element, these are based on the observed mean and
    # standard deviation in PCAWG tumors.
    mean.num.sigs.per.tumor = c(1.5, rep(3.9, 4)),
    sd.num.sigs.per.tumor   = c(0.5, rep(1.3, 4)),
    mean.log10.exp.per.sig  = c(3,   2.4, 2.1,  2.1, 2),
    sd.log10.exp.per.sig    = c(0.7, 0.6, 0.5,  0.4, 0.4)
  )
  return(retval)
}


#' Generate a suite of data with random signatures and exposures.
#'
#' @keywords internal
#'
#' Base parameters on SigProfiler attributions (exposures).

GenerateAllRandomSP <-
  function(top.level.dir  = "../SP.like.rand.exp.rand.sigs.2019.10.28",
           overwrite      = TRUE,
           verbose        = TRUE,
           num.replicates = 10,
           unlink         = FALSE) {
    log <- testthat::capture_messages(
      GenerateRandomExpRandomSigs200(
        parm                   = ParametersSPLike(),
        top.level.dir          = top.level.dir,

        # mean.log10.exp.per.sig = 2.97,
        # Based on mean(log10(sp.all.real.exposures[sp.all.real.exposures >= 1]))

        # sd.log10.exp.per.sig   = 0.7047,
        # Based on sd(log10(sp.all.real.exposures[sp.all.real.exposures >= 1]))

        seed                   = 1211,
        num.replicates         = num.replicates))
    if (!unlink) {
      cat(log, file = file.path(top.level.dir, "log.txt"))
    }
  }


# The next 2 functions provide the specifications for generating
# several suites of random data requested by Ludmil, 2019 08 25

Resummarize <- function(top.level.dir) {
  summary <- data.table::fread(file.path(top.level.dir, "summary.csv"))
  s1 <- as.data.frame(summary)
  fut.rownames <-
    paste("nsig", s1$total.num.sigs,
          "rep.num", s1$replicate.number, sep = ".")
  s1 <- s1[ , -ncol(s1)]
  rownames(s1) <- fut.rownames
  s2 <- split(s1, s1$total.num.sigs)
  add.means <- function(x) {
    xx <- colMeans(x)
    xxx <- rbind(x, as.list(xx))
    # xxxx <- as.character(signif(xxx, digits = 4))
    xxxx <-
      apply(xxx, MARGIN = 2, function(x) as.character(signif(x, digits = 3)))
    retv <- rbind(xxxx, "")
    return(retv)
  }

  s3 <- lapply(s2, add.means)
  s4 <- do.call("rbind", s3)
  write.csv(s4, file.path(top.level.dir, "summary2.csv"))
  return(s4)
}


GenerateRandomExpRandomSigs200 <- function(parm,
                                           top.level.dir,
                                           mean.log10.exp.per.sig,
                                           sd.log10.exp.per.sig,
                                           overwrite      = TRUE,
                                           verbose        = TRUE,
                                           num.replicates = 1,
                                           seed) {
  MustCreateDir(top.level.dir, overwrite)
  set.seed(seed, sample.kind = "Rejection")
  if (verbose) RNGMessages("GenerateRandomExpRandomSigs200")

  summary.file <- file.path(top.level.dir, "summary.csv")
  cat("total.num.sigs",
      "target.mean.num.sigs.per.tumor",
      "target.sd.num.sigs.per.tumor",
      "replicate.number",
      "actual.mean.num.sigs.per.tumor",
      "actual.sd.num.sigs.per.tumor",
      "mean.num.muts.per.tumor",
      "mean.log10.muts",
      "sd.num.muts.per.tumor",
      "sd.log10.muts",
      "median.num.muts.per.tumor",
      "mad.num.muts.per.tumor",
      "trail\n",
      sep = ",", file = summary.file)

  for (i in 1:nrow(parm)) {
    total.num.sigs          <- parm[i, "total.num.sigs"]
    mean.num.sigs.per.tumor <- parm[i, "mean.num.sigs.per.tumor"]
    sd.num.sigs.per.tumor   <- parm[i, "sd.num.sigs.per.tumor"]
    mean.log10.exp.per.sig  <- parm[i, "mean.log10.exp.per.sig"]
    sd.log10.exp.per.sig    <- parm[i, "sd.log10.exp.per.sig"]
    if (verbose) {
      message("\n", paste(rep("=", 40), collapse = ""))
      message("\nGenerateRandomExpRandomSigs200\n",
              "target distribution parameters:",
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

      stats.on.sim.data <-
        Create1CatRandomExpRandomSigs(
          num.syn.tumors          = 200,
          total.num.sigs          = total.num.sigs,
          mut.mean                = mean.log10.exp.per.sig,
          mut.sd                  = sd.log10.exp.per.sig,
          mean.num.sigs.per.tumor = mean.num.sigs.per.tumor,
          sd.num.sigs.per.tumor   = sd.num.sigs.per.tumor,
          dir.name                = dir.name,
          overwrite               = overwrite,
          verbose                 = verbose)

      cat(total.num.sigs, mean.num.sigs.per.tumor, sd.num.sigs.per.tumor,
          replicate.number, stats.on.sim.data,
          "\n", sep = ",", file = summary.file, append = TRUE)

    }
  }
}


#' Create a catalog of "random" synthetic spectra for 96-channel mutation types.
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
           verbose   = TRUE) {

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

    muts.per.tumor = colSums(exp)

    return(c(
      actual.sig.num.mean = attr(exp, "actual.sig.num.mean"),
      actual.sig.num.sd = attr(exp, "actual.sig.num.sd"),

      mean.num.muts.per.tumor = mean(muts.per.tumor),
      mean.log10.num = mean(log10(muts.per.tumor)),
      sd.num.muts.per.tumor = sd(muts.per.tumor),
      sd.log10.num = sd(log10(muts.per.tumor)),
      median.num.muts.per.tumor = stats::median(muts.per.tumor),
      mad.num.muts.per.tumor    = stats::mad(muts.per.tumor)

      ))
  }

