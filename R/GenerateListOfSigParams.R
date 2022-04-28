#' Generate a list of signature parameters for different cancer types from real exposure
#'
#' @inheritParams GenerateSyntheticTumors
#'
#' @note This function calls \code{GetSynSigParamsFromExposures}.
#'
#' @export
#'
#' @examples
#'
#' # Generate a list of signature parameters for Indel (ID) using negative binomial distribution
#' real.exposures.ID <- PCAWG7::exposure$PCAWG$ID
#' cancer.types <- PCAWG7::CancerTypes()[1:5]
#' sig.params <- SynSigGen::signature.params$ID
#' ID.sig.params <-
#'   GenerateListOfSigParams(real.exposures = real.exposures.ID,
#'                           cancer.types = cancer.types,
#'                           distribution = "neg.binom",
#'                           sig.params = sig.params
#'  )
GenerateListOfSigParams <- function(real.exposures, cancer.types, distribution = NULL,
                                    verbose = 0, sig.params = NULL)
{
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
  return(params)
}

