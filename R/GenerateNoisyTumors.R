#' Generate noisy tumors from available exposures
#'
#' @param seed A random seed to use.
#'
#' @param dir The directory in which to put the output; will be created if
#'   necessary.
#'
#' @param input.exposure A matrix of exposures.
#'
#' @param signatures A matrix of signatures.
#'
#' @param n.binom.size If non \code{NULL}, use negative binomial noise
#'     with this size parameter; see \code{\link[stats]{NegBinomial}}.
#'     If \code{NULL}, then use poisson distribution to do the resampling.
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
#' # Generate synthetic tumors for Indel (ID)
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
           n.binom.size = NULL, overwrite = TRUE) {
    # Set seed using R's default random number generator kind "Mersenne-Twister"
    set.seed(seed = seed, kind = "Mersenne-Twister")
    retval <- SynSigGen::AddNoise(input.exposure = input.exposure,
                                  signatures = signatures,
                                  n.binom.size = n.binom.size)

    if (overwrite == TRUE) {
      dir.create(path = dir, showWarnings = FALSE)
    } else {
      stop("\nDirectory ", dir, " exists\n")
    }

    ICAMSxtra::WriteExposure(exposure = retval$exposures,
                             file = file.path(dir,
                                              paste0("ground.truth.syn.exposures.noisy.n.binom.size.",
                                                     n.binom.size, ".csv")))
    ICAMS::WriteCatalog(catalog = ICAMS::as.catalog(retval$spectra),
                        file = file.path(dir,
                                         paste0("ground.truth.syn.catalog.noisy.n.binom.size.",
                                                n.binom.size, ".csv")))
    ICAMS::WriteCatalog(catalog = signatures,
                        file = file.path(dir, "ground.truth.syn.sigs.csv"))
    return(retval)
  }
