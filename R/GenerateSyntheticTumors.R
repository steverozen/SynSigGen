#' Generate synthetic tumors based on real exposures in one or more cancer types
#'
#' @param seed A random seed to use.
#'
#' @param dir The directory in which to put the output; will be created if
#'   necessary.
#'
#' @param cancer.types A vector of character strings denoting different cancer
#'   types. This function will search \code{real.exposures} for exposures from
#'   tumors matching these strings. See \code{PCAWG7::CancerTypes()} for example.
#'
#' @param samples.per.cancer.type Number of synthetic tumors to create for each
#'   cancer type. If it is \strong{one} number, then generate the \strong{same}
#'   number of synthetic tumors for each \code{cancer.types}. Or if it is a
#'   \strong{vector} of numbers, then generate synthetic tumors for each
#'   \code{cancer.type} accordingly to the number specified in the vector. The length
#'   and order of \code{samples.per.cancer.type} should match that in \code{cancer.types}.
#'
#' @param input.sigs A matrix of signatures.
#'
#' @param real.exposures A matrix of real exposures.
#'
#' @param distribution Probability distribution used to generate synthetic
#'   exposures due to active mutational signatures. Can be \code{neg.binom}
#'   which stands for negative binomial distribution. If \code{NULL} (Default),
#'   then this function uses log normal distribution with base 10.
#'
#' @param sample.prefix.name Prefix name to add to the synthetic tumors.
#'
#' @param tumor.marker.name Tumor marker name to add to the synthetic tumors.
#'   E.g. "MSI", "POLE".
#'
#' @param overwrite If TRUE, overwrite existing directories and files.
#'
#' @param verbose If > 0 cat various messages.
#'
#' @param sig.params Empirical signature parameters generated using real
#'   exposures irrespective of their cancer types. If there
#'   is only one tumor having a signature in a cancer type in \code{real.exposures},
#'   we cannot fit the \code{distribution} to only one data point. Instead, we
#'   will use the empirical parameter from \code{sig.params}.
#'   Users can use \code{SynSigGen:::GetSynSigParamsFromExposuresOld} to generate
#'   their own signature parameters. If \code{NULL}(default), this function uses the
#'   PCAWG7 empirical signature parameters. See \code{signature.params} for more details.
#'
#' @return A list of three elements that comprise the
#' synthetic data: \enumerate{
#'  \item \code{ground.truth.catalog}: Spectra catalog with rows denoting mutation
#'  types and columns denoting sample names.
#'  \item \code{ground.truth.signatures}: Signatures active
#'  in \code{ground.truth.catalog}.
#'  \item \code{ground.truth.exposures}: Exposures of \code{ground.truth.signatures}
#'  in \code{ground.truth.catalog}.
#' }
#'
#' @export
#'
#' @examples
#'
#' # Generate synthetic tumors for DBS78
#' input.sigs.DBS78 <- PCAWG7::signature$genome$DBS78
#' real.exposures.DBS78 <- PCAWG7::exposure$PCAWG$DBS78
#' cancer.types <- PCAWG7::CancerTypes()[1:5]
#' DBS78.synthetic.tumors <-
#'   GenerateSyntheticTumors(seed = 191906,
#'                           dir = file.path(tempdir(), "DBS78.synthetic.tumors"),
#'                           cancer.types = cancer.types,
#'                           samples.per.cancer.type = 30,
#'                           input.sigs = input.sigs.DBS78,
#'                           real.exposures = real.exposures.DBS78,
#'                           sample.prefix.name = "SP.Syn."
#'   )
#'
#' # Generate synthetic tumors for Indel (ID) using negative binomial distribution
#' input.sigs.ID <- PCAWG7::signature$genome$ID
#' real.exposures.ID <- PCAWG7::exposure$PCAWG$ID
#' cancer.types <- PCAWG7::CancerTypes()[1:5]
#' ID.synthetic.tumors <-
#'   GenerateSyntheticTumors(seed = 191906,
#'                           dir = file.path(tempdir(), "ID.synthetic.tumors"),
#'                           cancer.types = cancer.types,
#'                           samples.per.cancer.type = 30,
#'                           input.sigs = input.sigs.ID,
#'                           real.exposures = real.exposures.ID,
#'                           distribution = "neg.binom",
#'                           sample.prefix.name = "SP.Syn."
#'   )
#'
#' # Plot the synthetic catalog and exposures
#' ICAMS::PlotCatalogToPdf(catalog = DBS78.synthetic.tumors$ground.truth.catalog,
#'                         file = file.path(tempdir(), "DBS78.synthetic.catalog.pdf"))
#' ICAMSxtra::PlotExposureToPdf(exposure = DBS78.synthetic.tumors$ground.truth.exposures,
#'                              file = file.path(tempdir(), "DBS78.synthetic.exposures.pdf"),
#'                              cex.xaxis = 0.7)
GenerateSyntheticTumors <- function(seed,
                                    dir,
                                    cancer.types,
                                    samples.per.cancer.type,
                                    input.sigs,
                                    real.exposures,
                                    distribution = NULL,
                                    sample.prefix.name = "SP.Syn.",
                                    tumor.marker.name = NULL,
                                    overwrite       = TRUE,
                                    verbose = 0,
                                    sig.params = NULL)
{
  # Check whether the signatures in real.exposures are all available in input.sigs
  sigs.not.available <- setdiff(rownames(real.exposures), colnames(input.sigs))
  if (length(sigs.not.available) > 0) {
    stop("Some signatures in real.exposures are not available in input.sigs: ",
         paste(sigs.not.available, collapse = " "))
  }

  # Check whether there are samples in real.exposures that belong to the specified
  # cancer.types
  to.check <- lapply(cancer.types, FUN = function(x) {
    one.cancer.type <- x
    sample.exists <- any(grepl(pattern = one.cancer.type, x = colnames(real.exposures)))
    if(!sample.exists) {
      stop("Cannot find any sample in real.exposures that belong to cancer type: ",
           one.cancer.type,
           ". Please make sure the sample names in real.exposure have cancer type information.")
    }
  })

  # Check whether to generate different number of synthetic tumors for each cancer type
  if (length(samples.per.cancer.type) > 1) {
    if (length(cancer.types) != length(samples.per.cancer.type)) {
      stop("The number of values in samples.per.cancer.type should match the ",
           "number of elements in cancer.types")
    }
  } else {
    samples.per.cancer.type <- rep(samples.per.cancer.type, length(cancer.types))
  }

  # Get the mutation type of the synthetic data
  mutation.type <- GetMutationType(sig.name = colnames(input.sigs))

  suppressWarnings(set.seed(seed = seed, sample.kind = "Rounding"))

  # Getting empirical estimates of key parameters describing exposures due to signatures
  params <- lapply(cancer.types, FUN = function(x) {
    indices.one.type <- grep(pattern = x, x = colnames(real.exposures))
    exposures.one.type <- real.exposures[, indices.one.type, drop = FALSE]
    return(GetSynSigParamsFromExposures(exposures = exposures.one.type,
                                        distribution = distribution,
                                        verbose = verbose,
                                        cancer.type = x,
                                        sig.params = sig.params))
  })
  names(params) <- cancer.types

  # Generate synthetic exposures from real exposures
  synthetic.exposures <- lapply(1:length(cancer.types), FUN = function(x) {
    one.cancer.type <- cancer.types[x]
    params.one.type <- params[[one.cancer.type]]

    if (is.null(tumor.marker.name)) {
      sample.prefix.names <- paste0(sample.prefix.name, one.cancer.type, "::S")
    } else {
      sample.prefix.names <- paste0(sample.prefix.name, one.cancer.type, "::",
                                    tumor.marker.name, ".S")
    }

    # In case there is only one contributing signature to the cancer type
    if (ncol(params.one.type) == 1) {
      synthetic.exposures <-
        GenerateSyntheticExposures(sig.params = params.one.type,
                                   num.samples = samples.per.cancer.type[x],
                                   name = sample.prefix.names,
                                   sig.matrix = input.sigs,
                                   distribution = distribution)
      synthetic.exposures <- matrix(data = synthetic.exposures,
                                    ncol = samples.per.cancer.type[x],
                                    dimnames = list(colnames(params.one.type),
                                                    names(synthetic.exposures)))
      return(synthetic.exposures)

    } else {
      return(GenerateSyntheticExposures(sig.params = params.one.type,
                                        num.samples = samples.per.cancer.type[x],
                                        name = sample.prefix.names,
                                        sig.matrix = input.sigs,
                                        distribution = distribution))
    }

  })
  names(synthetic.exposures) <- cancer.types

  # Merge all exposure matrices in a list of matrices
  merged.exposures <- MergeExposures(synthetic.exposures)

  # Sort the signatures according to number ID
  merged.exposures.sorted.rowname <-
    merged.exposures[SortSigId(rownames(merged.exposures)), ]

  # Write parameters into files
  froot <- file.path(dir, "parameters")
  if (overwrite == TRUE) {
    dir.create(path = froot, showWarnings = FALSE, recursive = TRUE)
  } else {
    stop("\nDirectory ", froot, " exists\n")
  }

  lapply(cancer.types, FUN = function(x) {
    one.cancer.type <- x
    parms <- params[[one.cancer.type]]

    parm.file <- file.path(froot, paste0(sample.prefix.name, one.cancer.type,
                                         ".parms.", mutation.type, ".csv"))
    cat("# Original paramaters\n", file = parm.file)
    suppressWarnings( # Suppress warning on column names on append
      WriteSynSigParams(parms, parm.file, append = TRUE, col.names = NA))

    syn.exp <- synthetic.exposures[[one.cancer.type]]

    # Sanity check; we regenerate the parameters from the synthetic exposures.
    if (verbose > 0) {
      cat("\nRegenerating parameters from the synthetic exposures and",
          "compare with that from real exposures\n")
    }
    check.params <- GetSynSigParamsFromExposures(exposures = syn.exp,
                                                 distribution = distribution,
                                                 verbose = verbose,
                                                 cancer.type = one.cancer.type,
                                                 sig.params = sig.params)

    # check.params should be similar to parms
    cat("# Parameters derived from synthetic exposures\n",
        file = parm.file, append = TRUE)

    missing.sig.names <- setdiff(colnames(parms), colnames(check.params))
    if (length(missing.sig.names) > 0) {
      check.param2 <- matrix(NA, nrow = dim(parms)[1], ncol = dim(parms)[2])
      dimnames(check.param2) <- dimnames(parms)
      check.param2[ , colnames(check.params)] <- check.params
      check.params <- check.param2
    }

    suppressWarnings(
      WriteSynSigParams(check.params, parm.file, append = TRUE, col.names = NA))

    if (length(missing.sig.names) > 0) {
      cat("# Some signatures not represented in the synthetic data:\n",
          file = parm.file, append =  TRUE)
      cat("#", missing.sig.names, "\n", file = parm.file, append = TRUE)
    }

    cat("# Difference between original parameters and parameters",
        "derived from synthetic exposures\n",
        file = parm.file, append = TRUE)
    WriteSynSigParams(parms - check.params, parm.file,
                      append = TRUE, col.names = NA)
    if (is.null(distribution)) {
      cat("# The difference should be small\n",
          file = parm.file, append = TRUE)
    }
  }
  )

  catalog <- CreateAndWriteCatalog(sigs = input.sigs,
                                   exp = merged.exposures.sorted.rowname,
                                   my.dir = dir,
                                   overwrite = overwrite,
                                   extra.file.suffix = mutation.type)
  ground.truth.signatures <-
    input.sigs[, rownames(merged.exposures.sorted.rowname), drop = FALSE]
  return(list(ground.truth.catalog = catalog,
              ground.truth.signatures = ground.truth.signatures,
              ground.truth.exposures = merged.exposures.sorted.rowname))
}

#' @keywords internal
NumFromId<- function(s) {
  return(
    as.numeric(
      sub("[^0123456789]*(\\d+).*", "\\1", s, perl = TRUE)))
}

#' @keywords internal
SortSigId <- function(sig.id) {
  num <- NumFromId(sig.id)
  sig.id2 <- sig.id[order(num)]
  return(sig.id2)
}
