#' @title Write exposure matrix to a file
#'
#' @param exposure.matrix Matrix of exposures
#'
#' @param file File to which to write the exposure matrix (as a CSV file)
#'
#' @export
#'
#' @importFrom utils write.csv
#'
WriteExposure <- function(exposure.matrix, file) {
  old.digits <- getOption("digits")
  options(digits = 22)
  write.csv(exposure.matrix, file, row.names = TRUE)
  options(digits = old.digits)
}


#' @title Read an exposure matrix from a file
#'
#' @param file CSV file containing an exposure matrix
#'
#' @param check.names Passed to \code{utils::read.csv}.
#' IMPORTANT: If \code{TRUE} this will replace the double
#' colon in identifiers of the form <tumor_type>::<sample_id>
#' with two periods (i.e. <tumor_type>..<sample_id>.
#' If \code{check.names} is true, generate a warning
#' if double colons were present.
#'
#' @return Matrix of exposures
#'
#' @export
#'
#' @importFrom utils read.csv
ReadExposure <- function(file, check.names = TRUE) {
  if (check.names) {
    headers <- read.csv(file, nrow = 1, header = FALSE, stringsAsFactors = FALSE)
    double.colon <- grep("::", unlist(headers)[-1], fixed = TRUE)
    if (length(double.colon) > 0) {
      warning(":: in sample ID replaced by ..; suggest calling with check.names = FALSE")
    }
  }
  retval <- read.csv(file, row.names = 1, check.names = check.names)
  stopifnot(!duplicated(colnames(retval)))
  return(retval)
}


#' @title Read an exposure matrix from a Synapse file
#'
#' @param file CSV file containing an exposure matrix
#'
#' @return Matrix of exposures
#'
#' @export
#'
#' @importFrom utils read.csv read.table
ReadSynapseExposure <- function(file) {
  return(read.table(file, header = T, sep = "\t",
                    as.is = T, row.names =  1))
}
