#' Create synthetic exposure data.
#' @keywords internal

CreateRandomExposures <- function(num.exposures,
                                  mean.num.sigs.per.tumor,
                                  sd.num.sigs.per.tumor,
                                  total.num.sigs,
                                  per.sig.mean.and.sd,
                                  sample.name.prefix = "S",
                                  sigs,
                                  verbose) {

  buffer <- 100

  sig.nums <-
    CreateNumSigsPerTumor(
      num.exposures           = num.exposures + buffer,
      mean.num.sigs.per.tumor = mean.num.sigs.per.tumor,
      sd.num.sigs.per.tumor   = sd.num.sigs.per.tumor,
      total.num.sigs          = total.num.sigs)

  ss <- summary(sig.nums)
  actual.sd <- sd(sig.nums)
  actual.mean <- mean(sig.nums)
  if (verbose) {
      message("\n", paste(rep(".", 40), collapse = ""),
              "\nCreateRandomExposures, actual distribution:\n")
    for (nn in names(ss)) {
      message(nn, " = ", ss[nn])
    }
    message("sd = ", actual.sd, "\n")
  }

  exp <-
    sapply(sig.nums,
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

  attr(exp, "actual.sig.num.mean") <- actual.mean
  attr(exp, "actual.sig.num.sd")   <- actual.sd

  return(exp)
}


#' Create \code{num.exposures} signature counts from a normal distribution with \code{mean} and \code{sd}.
#'
#' Discard tumors with signature count = 0 or
#' signature count > \code{total.num.sigs}.
#'
#' @param num.exposures Number of exposures to create.
#'
#' @param mean.num.sigs.per.tumor Mean of distribution to draw from.
#'
#' @param sd Standard deviation of distribution to draw from.
#'
#' @param total.num.sigs Number of signatures in the "universe".
#'
#' @return A numeric vector, each element of which the number
#' of signatures in the corresponding tumor (which still remains
#' to be created.)
#'
#' @keywords internal

CreateNumSigsPerTumor <- function(num.exposures,
                                  mean.num.sigs.per.tumor,
                                  sd.num.sigs.per.tumor,
                                  total.num.sigs) {
  retval <- numeric(0)
  num.exposures.to.try <- num.exposures
  while (length(retval) < num.exposures) {
    num.exposures.to.try <- 3 * num.exposures.to.try
    retval <- round(rnorm(num.exposures.to.try,
                          mean = mean.num.sigs.per.tumor,
                          sd = sd.num.sigs.per.tumor))
    retval <- retval[retval > 0]
    retval <- retval[retval <= total.num.sigs]
  }
  return(retval[1:num.exposures])
}


#' Select \code{num.exp} signatures from the members of \code{sig.names} and create one column of an exposure matrix (as a vector).
#'
#' @param target.num.exp Number of signatures to which sample "was" exposed.
#'
#' @param all.sig.names Names of all possible signatures.
#'
#' @param target.sig.means Target means for creating mutation counts.
#'
#' @param target.sig.sds Target standard deviations for creating mutation counts.
#'
#' @return A subset of \code{sig.names} of size \code{num.exp}.
#'
#' @keywords internal

ExposureNums2Exposures <-
  function(target.num.exp, all.sig.names, target.sig.means, target.sig.sds) {
    stopifnot(names(target.sig.means) == names(target.sig.sds))
    stopifnot(all.sig.names == names(target.sig.means))
    all.len <- length(all.sig.names)
    index.to.use <- sample.int(all.len, size = target.num.exp, replace = FALSE)
    retval <- numeric(all.len) # all zeros
    retval[index.to.use] <-
      10^rnorm(target.num.exp,
               mean = target.sig.means[index.to.use],
               sd = target.sig.sds[index.to.use])
    names(retval) <- all.sig.names
    return(retval)
  }

