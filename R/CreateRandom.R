# This file contains functions to create "completely random"
# artificial signatures


#' This is the top-level function to create a set of spectra from random signatures.
#'
#' @param top.level.dir Directory in which to put all results. It will be
#' created if necessary.
#'
#' @param seed Use default for regression testing.
#'
#' @param regress.dir If not \code{NULL} compare the known results in
#'   this directory with the created results in \code{top.level.dir}.
#'
#' @param num.syn.tumors Total number of synthetic tumors to create. Use the
#' default for regression testing.
#'
#' @param overwrite If \code{TRUE} overwrite existing files and directories.
#'
#' @param unlink If \code{TRUE} unlink the created directory after the
#'   regression test.
#'
#' @param verbose If \code{TRUE} print a few informative messages.
#'
#' @export
#'
CreateRandomSyn <-
  function(
    top.level.dir,
    seed           = 1443196,
    regress.dir    = "data-raw/long.test.regression.data/syn.30.random.sigs/",
    num.syn.tumors = 1000,
    overwrite      = FALSE,
    unlink         = FALSE,
    verbose        = FALSE) {


    stopifnot(!missing(top.level.dir))

    suppressWarnings(RNGkind(sample.kind = "Rounding"))
    # For compatibility with R < 3.6.0
    set.seed(seed)

    CreateRandomSAAndSPSynCatalogs(top.level.dir  = top.level.dir,
                                   num.syn.tumors = num.syn.tumors,
                                   overwrite      = overwrite,
                                   verbose        = verbose)
    if (!is.null(regress.dir)) {
      return("ok" == NewDiff4SynDataSets(top.level.dir,
                                         regress.dir,
                                         unlink = unlink,
                                         verbose = TRUE))
     }
  }

# Test data for Ludmil and Mishu
create.1000.random.2019.08.30 <- function() {
  seeds <- sample(1000, size = 9)
  for (seed in seeds) {
    CreateRandomSyn(
      top.level.dir = paste0("../random.10000.", seed),
      seed = seed, overwrite = FALSE, regress.dir = NULL)
  }

}


#' Create a full SignatureAnalyzer / SigProfiler test data set for "random" artificial signatures.
#'
#' @param top.level.dir Path to top level of directory structure to be created.
#'
#' @param num.syn.tumors Number of synthetic tumors to create.
#'
#' @param overwrite If \code{TRUE}, overwrite existing directories / files.
#'
#' @keywords internal

CreateRandomSAAndSPSynCatalogs <-
  function(top.level.dir, num.syn.tumors, overwrite = FALSE, verbose = FALSE) {

    COMPOSITE.features <- c(ICAMS::catalog.row.order[["SBS1536"]],
                            ICAMS::catalog.row.order[["DBS78"]],
                            ICAMS::catalog.row.order[["ID"]])
    stopifnot(length(COMPOSITE.features) == 1697)

    MustCreateDir(top.level.dir, overwrite)

    # The following are for choosing the mean number of mutations due to each
    # synthetic signature.
    sa.mut.mean <- 2.349
    # Based on mean(log10(sa.all.real.exposures[sa.all.real.exposures >= 1]))

    sa.mut.sd   <- 0.6641
    # Based on sd(log10(sa.all.real.exposures[sa.all.real.exposures >= 1]))

    sp.mut.mean <- 2.97
    # Based on mean(log10(sp.all.real.exposures[sp.all.real.exposures >= 1]))

    sp.mut.sd   <- 0.7047
    # Based on sd(log10(sp.all.real.exposures[sp.all.real.exposures >= 1]))

    # An alternative would be:
    # per.sig.mean <- apply(sa.all.real.exposures, 1, function(x) mean(log10(x[ x > 1])
    # sa.mut.mean  <- mean(per.sig.mean, na.rm = TRUE)
    # sa.mut.sd    <- sd(per.sig.mean, na.rm = TRUE)

    # These are the mean and SD of the number of signatures per tumor.
    # E.g. sa.all.real.exposures > 0 return TRUE or FALSE, which
    # are coerced to 1 and 0 by mean and sd.
    sa.num.sigs.mean <- 15.525 # mean(colSums(sa.all.real.exposures > 0))
    sa.num.sigs.sd   <-  6.172 # sd(colSums(sa.all.real.exposures > 0))
    sp.num.sigs.mean <-  3.947 # mean(colSums(sp.all.real.exposures > 0))
    sp.num.sigs.sd   <-  1.331 # sd(colSums(sp.all.real.exposures > 0))

    num.sigs.to.create <- 30 # Also tried, 60 (ncol(SynSig::sa.96.sigs));
    # this is too many.

    CreateOnePairOfRandomCatalogs(
      num.syn.tumors     = num.syn.tumors,
      total.num.sigs     = num.sigs.to.create,
      mut.mean           = sa.mut.mean,
      mut.sd             = sa.mut.sd,
      num.sigs.mean      = sa.num.sigs.mean,
      num.sigs.sd        = sa.num.sigs.sd,
      sig.name.prefix    = "SARandSig",
      sample.name.prefix = "SARandSample",
      composite.dir.name = file.path(top.level.dir, "sa.sa.COMPOSITE"),
      x96.dir.name       = file.path(top.level.dir, "sa.sa.96"),
      COMPOSITE.features = COMPOSITE.features,
      overwrite          = overwrite,
      verbose            = verbose)

    CreateOnePairOfRandomCatalogs(
      num.syn.tumors     = num.syn.tumors,
      total.num.sigs     = num.sigs.to.create,
      mut.mean           = sp.mut.mean,
      mut.sd             = sp.mut.sd,
      num.sigs.mean      = sp.num.sigs.mean,
      num.sigs.sd        = sp.num.sigs.sd,
      sig.name.prefix    = "SPRandSig",
      sample.name.prefix = "SPRandSample",
      composite.dir.name = file.path(top.level.dir, "sp.sa.COMPOSITE"),
      x96.dir.name       = file.path(top.level.dir, "sp.sp"),
      COMPOSITE.features = COMPOSITE.features,
      overwrite          = overwrite,
      verbose            = verbose)

  }


#' Create one "random" artificial signature profile.
#'
#' @param row.names One of the \code{\link{ICAMS}} package variable such as
#'  \code{catalog.row.order[["SBS96"]]}.
#'
#' @return A single column matrix with \code{rownames} \code{row.headers} and
#'   \code{colnames} \code{"RandSig"}.
#'
#' @importFrom stats runif
#'
#' @keywords internal

CreateOneRandomMutSigProfile <- function(row.names) {
  stopifnot(!is.null(row.names))

  retval <- matrix(10^runif(length(row.names)), ncol = 1)
  # retval <- matrix(10^rnorm(length(row.names)), ncol = 1) # Too spiky
  retval <- retval / sum(retval)
  rownames(retval) <- row.names
  colnames(retval) <- "RandSig"
  return(retval)
 }

#' Create a matrix of "random" signature profiles.
#'
#' @param row.headers One of the \code{\link{ICAMS}} package variable such as
#'  \code{catalog.row.order[["SBS96"]]}.
#'
#' @param num.signatures Number of signatures to create.
#'
#' @param sig.name.prefix The signatures will be named \code{<sig.name.prefix>1},
#' \code{<sig.name.prefix>2}, etc.
#'
#' @return A \code{num.signatures}-column
#'   matrix with rownames \code{row.headers}.
#'
#' @keywords internal

CreateRandomMutSigProfiles <-
  function(row.headers, num.signatures, sig.name.prefix) {

  stopifnot(!is.null(row.headers))

  retval <- lapply(1:num.signatures,
                function(x) CreateOneRandomMutSigProfile(row.headers))
  retval <- as.matrix(data.frame(retval))
  colnames(retval) <- paste0(sig.name.prefix, 1:num.signatures)
  return(retval)
}

#' Create means and standard deviations of log10 mutation counts for synthetic signatures.
#'
#' @param num.sigs Number of signatures for which means and standard deviations
#' are needed.
#'
#' @param target.mut.mean Target number of mean of log10 of number
#'  of mutations per tumor due to one signature.
#'
#' @param target.mut.sd Target of standard deviation of  the the log10 of the
#' number of mutations per tumor due to one signature.
#'
#' @param sig.names Names for the output vectors.
#'
#' @keywords internal
#'
#' @return A list with the needed synthetic means
#' and standard deviations.

CreateMeanAndStdevForSigs <-
  function(num.sigs, target.mut.mean, target.mut.sd, sig.names) {
  syn.mean <- rnorm(num.sigs, mean = target.mut.mean, sd = target.mut.sd)
  names(syn.mean) <- sig.names
  syn.sd  <- syn.mean * target.mut.sd / target.mut.mean
  names(syn.sd) <- sig.names
  return(list(syn.mean = syn.mean, syn.sd = syn.sd))
}

#' Create a pair of "random" synthetic catalogs, one for 96-channel and one for COMPOSITE features.
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
#' @param num.sigs.mean Mean number of signatures contributing to each tumor.
#'
#' @param num.sigs.sd Standard deviation the number of signatures
#' contribution to each tumor.
#'
#' @param sig.name.prefix String to put in front of an integer (as
#' string) to form an identifier for a synthetic signature.
#'
#' @param sample.name.prefix String to put in front of an integer (as
#' string) to form an identifier for a synthetic sample (tumor).
#'
#' @param composite.dir.name string indicating the name of the COMPOSITE
#' subdirectory; probably one of \code{"sa.sa.COMPOSITE"} or
#' \code{"sp.sa.COMPOSITE"}.
#'
#' @param x96.dir.name A string indicating the name of the 96-channel
#' subdirectory; probably one of \code{"sa.sa.COMPOSITE"} or
#' \code{"sp.sa.COMPOSITE"}.
#'
#' @param COMPOSITE.features Character vector containing
#' rownames for a COMPOSITE signature or catalog.
#'
#' @param overwrite If \code{TRUE} overwrite existing directories / files.
#'
#' @keywords internal

CreateOnePairOfRandomCatalogs <-
  function(num.syn.tumors,
           total.num.sigs,
           mut.mean,
           mut.sd,
           num.sigs.mean,
           num.sigs.sd,
           sig.name.prefix,
           sample.name.prefix,
           composite.dir.name,
           x96.dir.name,
           COMPOSITE.features,
           overwrite = FALSE,
           verbose = TRUE) {

    syn.96.sigs <-
      CreateRandomMutSigProfiles(
        ICAMS::catalog.row.order[["SBS96"]], total.num.sigs, sig.name.prefix)

    syn.COMPOSITE.sigs <-
      CreateRandomMutSigProfiles(
        COMPOSITE.features, total.num.sigs, sig.name.prefix)

    sig.info <- CreateMeanAndStdevForSigs(
      total.num.sigs, mut.mean, mut.sd, colnames(syn.96.sigs))

    exp <- CreateRandomExposures(
        num.exposures = num.syn.tumors,
        mean.num.sigs.per.tumor = num.sigs.mean,
        sd.num.sigs.per.tumor = num.sigs.sd,
        total.num.sigs = total.num.sigs,
        per.sig.mean.and.sd = sig.info,
        sample.name.prefix = sample.name.prefix,
        sigs = syn.COMPOSITE.sigs,
        # Todo(Steve): The reason for using COMPOSITE is subtle -- need to review and document.
        verbose = verbose)

    NewCreateAndWriteCatalog(
      sigs      = syn.COMPOSITE.sigs,
      exp       = exp,
      dir       = composite.dir.name,
      overwrite = overwrite)

    NewCreateAndWriteCatalog(
      sigs      = syn.96.sigs,
      exp       = exp,
      dir       = x96.dir.name,
      overwrite = overwrite)


  }
