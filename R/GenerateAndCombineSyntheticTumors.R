GenerateAndCombineSyntheticTumors <- function(seed,
                                              dir,
                                              list.of.cancer.types,
                                              list.of.samples.per.cancer.type,
                                              input.sigs,
                                              list.of.real.exposures,
                                              distribution = NULL,
                                              sample.prefix.name = "SP.Syn.",
                                              overwrite       = TRUE,
                                              verbose = 0)
{
  # Check whether the length of several lists are the same
  length1 <- length(list.of.cancer.types)
  length2 <- length(list.of.samples.per.cancer.type)
  length3 <- length(list.of.real.exposures)

  if (length1 != length2) {
    stop("The number of elements in list.of.cancer.types should be the same ",
         "as that in list.of.samples.per.cancer.type")
  }

  if (length1 != length3) {
    stop("The number of elements in list.of.cancer.types should be the same ",
         "as that in list.of.real.exposures")
  }

  if (length2 != length3) {
    stop("The number of elements in list.of.samples.per.cancer.type should ",
         "be the same as that in list.of.real.exposures")
  }

  syn.exposures <- lapply(1:length1, FUN = function(x) {
    if (verbose > 0) {
      cat("Generating synthetic exposures for set ", x)
    }
    tmpdir <- tempfile()
    dir.create(tmpdir)
    retval <-
      GenerateSyntheticTumors(seed = seed,
                              dir = tmpdir,
                              cancer.types = list.of.cancer.types[[x]],
                              samples.per.cancer.type = list.of.samples.per.cancer.type[[x]],
                              input.sigs = input.sigs,
                              real.exposures = list.of.real.exposures[[x]],
                              distribution = distribution,
                              sample.prefix.name = sample.prefix.name,
                              overwrite = overwrite,
                              verbose = verbose)
    unlink(x = tmpdir, recursive = TRUE)
    return(retval$ground.truth.exposures)
  })



}

MergeExposuresByCancerTypes <- function(list.of.exposures) {
  exposures.by.cancer.types <- lapply(list.of.exposures, FUN = function(x) {
    one.exposure <- x
    return(PCAWG7::SplitPCAWGMatrixByTumorType(one.exposure))
  })

  lengths <- sapply(exposures.by.cancer.types, FUN = length)
  max.length.index <- which.max(lengths)
  max.num.of.cancer.types <- names(exposures.by.cancer.types[[max.length.index]])

  retval <- lapply(max.num.of.cancer.types, FUN = function(x) {
    one.cancer.type <- x

  })
}
