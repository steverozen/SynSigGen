#' @title Using parameters given to generate exposures for synthetic tumors
#'
#' @param tumor Signature presence matrix or exposure matrix for a tumor.
#' It has only one row, and K (# of signatures) columns.
#' Value in each column is the presence flag for a mutational signature:
#' the value can be non-zero(signature is present) or 0(absent).
#' The name of each column should be the name of a signature.
#'
#' @param sig.interest Names of mutational signatures you want to use to
#' generate exposures. It can be all, or part of signatures in colnames(tumor).
#'
#' @param distribution Probability distribution used to generate exposures for
#'   synthetic tumors. Can be either \code{neg.binom} (negative binomial) or
#'   \code{log.norm} (log normal).
#'
#' @param burden.per.sig Mean mutation burden of the
#' counts of mutations per megabase.
#' It has one row, and K columns. Each column name refers to a mutational signature.
#'
#' @param sd.per.sig standard deviation of mutation burden.
#' It has one row, and K columns. Each column name refers to a mutational signature.
#'
#' @details Determine the intensity of each
#' mutational signature in a tumor, returning the number of mutations
#' using the mean mutation burden per signature and the std dev
#'
#' @importFrom stats rnorm
#'
#' @keywords internal
GenerateSynExposureOneSampleNew <-
  function(tumor,
           sig.interest,
           distribution,
           burden.per.sig,
           sd.per.sig
  ) {

    if (!distribution %in% c("neg.binom", "log.norm")) {
      stop("Can only use neg.binom or log.norm for argument distribution")
    }

    ## starts with individual tumors, only generate exposures for signatures
    ## with a flag does not equal to 0.
    active.sigs <- base::which(tumor != 0)

    for (sigs in active.sigs) {
      stdev <- sd.per.sig[,sigs]
      burden <- burden.per.sig[,sigs]

      if (distribution == "neg.binom") {
        tumor[sigs] <- stats::rnbinom(n = 1, mu = burden, size = )

      } else {
        ## if std dev is too big, >= 3, max = 3
        ### consider handling this different. the worry is that the variation
        ##  is too large, the sampled mutation burden will be very high,
        ### which will have a mutation burden that is not biologically possible
        if (stdev >= 3) {
          cat("Very large stdev", stdev, "\n")
          stdev = 3
        }

        ## mutational intensity follows a log normal distibution
        ## use the normal distribution with log-ed values instead
        tumor[sigs] <- 10^(rnorm(1, sd = stdev, mean = burden))
      }
    }

    tumor <- as.matrix(tumor)
    names(tumor) <- sig.interest
    return(tumor)

  }

#' @title Create synthetic exposures based given parameters
#'
#' @return A matrix with the rows being each signature and the columns being
#' generated samples. Each entry is the count of mutations due to one
#' signature in one sample.
#'
#' @param sig.params Parameters from \code{\link{GetSynSigParamsFromExposures}}
#'   or another source. Should be
#'   a matrix or data frame with one column for
#'   each signature and the following rows:
#' \describe{
#' \item{prob}{The proportion of tumors with the signature.}
#' \item{mean}{The mean(log_10(number of mutations)).}
#' \item{stdev}{The stdev(log_10(number of mutations)).}
#' }
#'   The rownames need to be the column names of a signature
#'   catalog.
#'
#' @param num.samples Number of samples to generate
#'
#' @param name Prefix for sample identifiers in the simulated dataset
#'
#' @export

GenerateSyntheticExposuresNew <-
  function(sig.params,
           num.samples = 10,
           name = 'synthetic') {

    sigs <- colnames(sig.params)
    stopifnot(!is.null(sigs))
    prev.present <- unlist(sig.params['prob', ]) # Note, get a vector
    sig.burden <- sig.params['mean', , drop = FALSE]
    sig.sd <- sig.params['stdev', , drop = FALSE]

    sig.present <- present.sigs(num.samples, prev.present)

    colnames(sig.present) <- paste(name, seq(1, num.samples), sep = '.')

    # Create a synthetic exposures for each column (sample)
    # in sig.present.
    retval <-
      apply(sig.present,
            2,
            GenerateSynExposureOneSample,
            sigs,
            sig.burden, ## burden is in mutation per megabase
            sig.sd)
    return(retval)
  }



#' Create a test data set based on >= 1 tumor types.
#'
#' @param top.level.dir Path to top level of directory structure to be created.
#'
#' @param cancer.type.strings Search the PCAWG data for tumors matching
#' these strings. Each string should identify one tumor type, for
#' some definition of tumor type. Probably the tumors in each type
#' should be non-overlapping, but the code does not enforce this and
#' does not care.
#'
#' @param num.syn.tumors Number of synthetic tumors to create
#' for each cancer type.
#'
#' @param overwrite If TRUE, overwrite existing directories / files.
#'
#' @param sa.exp SignatureAnalyzer exposures from which to select cancer types
#'        specified by \code{cancer.type.strings}. In the
#'        column names of \code{sa.exp} the cancer type string
#'        should be separated from the sample identifier by two colons
#'        (::).
#'
#' @param sp.exp SigProfiler exposures from which to select cancer types
#'        specified by \code{cancer.type.strings}. In the
#'        column names of \code{sp.exp} the cancer type string
#'        should be separated from the sample identifier by two colons
#'        (::).
#
#' @param verbose If > 0, cat various messages.
#'
#' @param bladder.regress.hack For use by \code{\link{BladderSkin1000}}.
#'     Forces use of non-hyper-mutated exposures for bladder-TCC even if
#'     \code{sa.exp} and \code{sp.exp} include hyper-mutated exposures.
#'
#' @export

CreateMixedTumorTypeSyntheticDataNew <- function(top.level.dir,
                                              cancer.type.strings,
                                              num.syn.tumors,
                                              overwrite = FALSE,
                                              sa.exp = sa.all.real.exposures,
                                              sp.exp = sp.all.real.exposures,
                                              verbose = FALSE,
                                              bladder.regress.hack = FALSE) {

  odigits <- getOption("digits")
  options(digits = 9)
  on.exit(options(digits = odigits))

  info.list <-
    lapply(cancer.type.strings,
           function(ca.type.str) {
             if (verbose) message("\n\nProcessing", ca.type.str, "\n\n\n")
             local.sa.exp <- sa.exp
             local.sp.exp <- sp.exp
             if (bladder.regress.hack && ca.type.str == "Bladder-TCC") {
               local.sa.exp <- SynSigGen::sa.no.hyper.real.exposures
               local.sp.exp <- SynSigGen::sp.no.hyper.real.exposures
               message("bladder.regress.hack deployed for ", ca.type.str)
             }
             retval <-
               SAAndSPSynDataOneCAType(
                 sa.real.exp    = local.sa.exp,
                 sp.real.exp    = local.sp.exp,
                 ca.type        = ca.type.str,
                 num.syn.tumors = num.syn.tumors,
                 file.prefix    = ca.type.str,
                 top.level.dir  = top.level.dir)
             return(retval)
           })

  # info.list has both SignatureAnalyzer and
  # SigProfiler exposures.

  sa.exposures <- lapply(info.list, function(x) x$sa.syn.exp)
  sa.exp <- MergeExposures(sa.exposures)
  sa.exp <- sa.exp[rowSums(sa.exp) > 0.5, ]
  if (verbose) message("Dimension sa.exp", dim(sa.exp))

  sp.exposures <- lapply(info.list, function(x) x$sp.syn.exp)
  sp.exp <- MergeExposures(sp.exposures)
  sp.exp <- sp.exp[rowSums(sp.exp) > 0.5, ]
  if (verbose) message("Dimension sp.exp", dim(sp.exp))


  # We will need the exposures later when evaluating the attributed signatures
  WriteExposure(sa.exp, file.path(top.level.dir, "sa.exposure.csv"))
  WriteExposure(sp.exp, file.path(top.level.dir, "sp.exposure.csv"))

  # Create catalogs of synthetic mutational spectra
  # based on SignatureAnalyzer attributions

  CreateAndWriteCatalog(
    sa.COMPOSITE.sigs,
    sa.exp,
    dir = NULL, # "sa.sa.COMPOSITE",
    WriteCatCOMPOSITE,
    overwrite = overwrite,
    my.dir = file.path(top.level.dir, "sa.sa.COMPOSITE"))

  CreateAndWriteCatalog(
    sa.96.sigs,
    sa.exp,
    dir = NULL, # "sa.sa.96",
    WriteCatalog,
    overwrite = overwrite,
    my.dir = file.path(top.level.dir, "sa.sa.96"))

  # Create catalogs of synthetic mutational spectra
  # based on SigProfiler attributions

  # First we need the matching between SigProfiler and
  # SignatureAnalyzers signatures.

  sp.sa.map.info <-
    MapSPToSASignatureNamesInExposure(sp.exp)

  CreateAndWriteCatalog(
    sa.COMPOSITE.sigs,
    sp.sa.map.info$exp2,
    dir = NULL, # "sp.sa.COMPOSITE",
    WriteCatCOMPOSITE,
    overwrite = overwrite,
    my.dir = file.path(top.level.dir, "sp.sa.COMPOSITE"))


  if (verbose) print(sp.sa.map.info$sp.to.sa.sig.match)

  CreateAndWriteCatalog(
    sp.sigs,
    sp.exp,
    dir = NULL, # "sp.sp",
    WriteCatalog,
    overwrite = overwrite,
    my.dir = file.path(top.level.dir, "sp.sp"))

  # AddAllScripts(maxK = 50)

  invisible(list(info.list = info.list,
                 sp.sa.map.info =  sp.sa.map.info))

}


#' @title Empirical estimates of key parameters describing exposures due to signatures.
#'
#' @param exposures A matrix in which each column is a sample and each row is a mutation
#'         signature, with each element being the "exposure",
#'         i.e. mutation count attributed to a
#'         (sample, signature) pair.
#'
#' @param verbose If > 0 cat various messages.
#'
#' @return A data frame with one column for
#' each of a subset of the input signatures
#' and the following rows
#' \enumerate{
#' \item the proportion of tumors with the signature
#' \item mean(log_10(mutations.per.Mb))
#' \item stdev(log_10(mutations.per.Mb))
#' }
#' Signatures not present in
#'  \code{exposures} or present only in a single tumor in
#'  \code{exposures} are removed.
#'
#' @export

GetParamsFromExposuresNew <- function(exposures, verbose = 0) {
  stopifnot(ncol(exposures) > 0)

  integer.counts <- round(exposures, digits = 0)
  integer.counts <- integer.counts[rowSums(integer.counts) > 0 , ]
  ret1 <- apply(X      = integer.counts,
                MARGIN = 1,
                FUN    = SynSigParamsOneSignature)

  # Some standard deviations can be NA (if there is only one tumor
  # with mutations for that signature). We pretend we did not see
  # these signatures. TODO(Steve): impute from similar signatures.
  if (any(is.na(ret1['stdev', ]))) {
    if (verbose > 0) {
      cat("\nWarning, some signatures present in only one sample, dropping:\n")
      cat(colnames(ret1)[is.na(ret1['stdev', ])], "\n")
    }
  }
  retval <- ret1[,!is.na(ret1['stdev',]) , drop = FALSE]
  if (ncol(retval) == 0) {
    stop("No signatures with usable parameters (> 1 sample with exposure)")
  }
  return(retval)
}


GenerateSyntheticTumorsNew <- function() {

}
