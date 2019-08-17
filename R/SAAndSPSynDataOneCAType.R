#' Generate parallel synthetic exposures from real SA and SP exposures and signatures
#'
#' @export
#'
#' @param sa.real.exp Exposure matrix from SignatureAnalyzer.
#'
#' @param sp.real.exp Exposure matrix from SigProfiler.
#'
#' @param ca.type The type the cancer, which is used in sample identifiers,
#'   which SigProfiler expects.
#'
#' @param num.syn.tumors Number of synthetic tumors to generate.
#'
#' @param file.prefix To explain later.
#'
#' @param top.level.dir Specifies the location to generate files.
#'
#' @return A list with the following elements:
#'
#' \enumerate{
#'
#'  \item \code{sa.parms} The parameters computed from \code{sa.real.exp}.
#'        This a matrix with a column for each signature
#'        and 3 rows:
#'
#'  \enumerate{
#'
#'    \item The proportion of tumors with given signature (in
#'          \code{sa.real.exp}).
#'
#'    \item The mean of the log10 of the number of mutations for a
#'          given siganture.
#'
#'    \item The standard deviation of log10 of the number of mutations
#'          for a given signature.
#'  }
#'
#'  \item \code{sa.syn.exp} The synthetic exposures computed from
#'        \code{sa.parms}.
#'
#'  \item \code{sp.parms} The parameters computed from \code{sp.real.exp}, with
#'         rows analogous to the rows in \code{sa.parms}.
#'
#'  \item \code{sp.syn.exp} The synthetic exposures computed from
#'        \code{sp.parms}.
#'
#'  }
#'
#'  @details Creates a bunch of files in location
#'  governed by \code{top.level.dir}. The main rationale for packaging this
#'  as one function is to ensure that some conventions regarding file
#'  naming are followed.
#'
#'  This function does \strong{not} create the synthetic
#'  mutational spectra catalogs but \strong{does} generate the
#'  synthetic exposures.

SAAndSPSynDataOneCAType <-
  function(sa.real.exp,
           sp.real.exp,
           ca.type,
           num.syn.tumors,
           file.prefix,
           top.level.dir = NULL) {


    sp.real.exp <- sp.real.exp[ , colnames(sa.real.exp)]

    stopifnot(colnames(sp.real.exp) == colnames(sa.real.exp))


    sa.real.exp <- GetExpForOneCancerType(ca.type, sa.real.exp)
    sp.real.exp <- GetExpForOneCancerType(ca.type, sp.real.exp)
    # ca.type <- paste0(ca.type, "::")
    # samples.to.use <-
    #  grep(ca.type, colnames(sa.real.exp), fixed = TRUE)
    # sa.real.exp <- sa.real.exp[ , samples.to.use]
    # sp.real.exp <- sp.real.exp[ , samples.to.use]


    # stopifnot(sa.real.exp == sa2)
    # stopifnot(sp.real.exp == sp2)

    # Make sure the 2 exposure data sets have the same colummns
    # in the same order.
    stopifnot(colnames(sa.real.exp) == colnames(sp.real.exp))

    sa.info <-
      GenerateSynFromReal(
        sa.real.exp,
        num.syn.tumors,
        file.prefix = paste0("sa.", file.prefix),
        sample.id.prefix = paste0("SA.Syn.", ca.type, "::S"),
        top.level.dir = top.level.dir)
    sp.info <-
      GenerateSynFromReal(
        sp.real.exp,
        num.syn.tumors,
        file.prefix = paste0("sp.", file.prefix),
        sample.id.prefix = paste0("SP.Syn.", ca.type, "::S"),
        top.level.dir = top.level.dir)
    return(list(
      sa.parms  = sa.info$parms,
      sa.syn.exp = sa.info$syn.exp,
      sp.parms   = sp.info$parms,
      sp.syn.exp = sp.info$syn.exp))
  }

GetExpForOneCancerType <- function(ca.type, exp) {
  ca.type <- paste0(ca.type, "::")
  samples.to.use <- grep(ca.type, colnames(exp), fixed = TRUE)
  return(exp[ , samples.to.use, drop = FALSE])
}
