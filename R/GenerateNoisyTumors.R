#' Generate noisy tumors from available exposures
#'
#' @inheritParams AddNoise
#'
#' @param seed A random seed to use.
#'
#' @param dir The directory in which to put the output; will be created if
#'   necessary.
#'
#' @param overwrite If TRUE, overwrite existing directories and files.
#'
#' @return A list with the elements \describe{
#' \item{expsoures}{The numbers of mutations due to each signature
#'    after adding noise}
#' \item{spectra}{The spectra based on the noisy signature exposures.}
#' }
#'
#' @export
#'
#' @examples
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
#' # Add noise to the exposures
#' ID.noisy.tumors <-
#'   GenerateNoisyTumors(seed = 892513,
#'                       dir = file.path(tempdir(), "ID.noisy.tumors"),
#'                       input.exposure = ID.synthetic.tumors$ground.truth.exposures,
#'                       signatures = ID.synthetic.tumors$ground.truth.signatures,
#'                       n.binom.size = 1)
#'
#' # Plot the synthetic and noisy catalog and exposures
#' ICAMS::PlotCatalogToPdf(catalog = ID.synthetic.tumors$ground.truth.catalog,
#'                         file = file.path(tempdir(), "ID.synthetic.catalog.pdf"))
#' ICAMSxtra::PlotExposureToPdf(exposure = ID.synthetic.tumors$ground.truth.exposures,
#'                              file = file.path(tempdir(), "ID.synthetic.exposures.pdf"),
#'                              cex.xaxis = 0.7)
#' ICAMS::PlotCatalogToPdf(catalog = ID.noisy.tumors$spectra,
#'                         file = file.path(tempdir(), "ID.noisy.catalog.pdf"))
#' ICAMSxtra::PlotExposureToPdf(exposure = ID.noisy.tumors$exposures,
#'                              file = file.path(tempdir(), "ID.noisy.exposures.pdf"),
#'                              cex.xaxis = 0.7)
#'
GenerateNoisyTumors <-
  function(seed, dir, input.exposure, signatures,
           n.binom.size = NULL, cp.factor = NULL, overwrite = TRUE) {
    # Set seed using R's default random number generator kind "Mersenne-Twister"
    set.seed(seed = seed, kind = "Mersenne-Twister")
    retval <- SynSigGen::AddNoise(input.exposure = input.exposure,
                                  signatures = signatures,
                                  n.binom.size = n.binom.size,
                                  cp.factor = cp.factor)

    if (overwrite == TRUE) {
      dir.create(path = dir, showWarnings = FALSE)
    } else {
      stop("\nDirectory ", dir, " exists\n")
    }

    # Get the mutation type of the noisy data
    mutation.type <- GetMutationType(sig.name = colnames(signatures))

    if (!is.null(cp.factor)) {
      exposure.file.suffix <- ""
      catalog.file.suffix <-
        paste0(".dirichlet.multinom.noise.cp.factor.", cp.factor)
    } else if (!is.null(n.binom.size)) {
      exposure.file.suffix <- paste0(".noisy.neg.binom.size.", n.binom.size)
      catalog.file.suffix <- paste0(".noisy.neg.binom.size.", n.binom.size)
    } else {
      exposure.file.suffix <- ".poisson.noise"
      catalog.file.suffix <- ".poisson.noise"
    }

    ICAMSxtra::WriteExposure(exposure = retval$exposures,
                             file = file.path(dir,
                                              paste0("ground.truth.syn.exposures",
                                                     exposure.file.suffix,
                                                     ".", mutation.type, ".csv")))
    ICAMS::WriteCatalog(catalog = ICAMS::as.catalog(retval$spectra),
                        file = file.path(dir,
                                         paste0("ground.truth.syn.catalog",
                                                catalog.file.suffix,
                                                ".", mutation.type, ".csv")))
    ICAMS::WriteCatalog(catalog = signatures,
                        file = file.path(dir,
                                         paste0("ground.truth.syn.sigs",
                                                ".", mutation.type, ".csv")))
    return(retval)
  }
