#' Add an R script to run SignatureAnalyzer a particular subdirectory.
#'
#' @param maxK The \code{maxK} argument for SignatureAnalyzers.
#'
#' @param slice Which subdirectory to put the script into.
#'
#' @param dir.name The name of the subdirectory.
#'
#' @keywords internal

AddScript <- function(maxK, slice, dir.name) {

  lines <- c(
    "",
    "",
    "",
    "library(SynSig)",
    "library(ICAMS)",
    "cat(\"\\n\\nRunning, maxK.for.SA is\", maxK.for.SA, \"\\n\\n\")",
    "RNGkind(kind = \"L'Ecuyer-CMRG\")",
    "set.seed(888)",
    "",
    "reval <- SignatureAnalyzer4MatchedCatalogs(",
    "  num.runs = 20,",
    "  signatureanalyzer.code.dir = \"/home/gmssgr/bin/SignatureAnalzyer.052418/\",",
    "  dir.root = \"..\",",
    "",
    "  overwrite = FALSE,",
    "  maxK = maxK.for.SA,",
    "  mc.cores = 20",
    "  )"
  )

  out.script.name <- paste0(slice, ".run.SA.R")
  lines[1]  <-
    paste0("# Put this file in <top.level.dir>/", dir.name,
           " and run Rscript ", out.script.name)
  lines[2]  <- paste0("maxK.for.SA <- ", maxK)
  lines[14] <- paste0("  slice = ", slice, ",")
  out.name  <- paste0(dir.name, "/", out.script.name)
  writeLines(lines, con = out.name)
}


#' Create scripts to run SignatureAnalyzer in all subdirectories of \code{top.level.dir}.
#'
#' WARNING, not tested recently as of 2020 05 12
#'
#' @param maxK The \code{maxK} argument for SignatureAnalyzer.
#'
#' @param top.level.dir Add the scripts to sub-directories
#'    \code{sa.sa.96}, \code{sp.sp}, etc. of this
#'    directory.
#'
#' @keywords internal

AddAllScripts <- function(maxK = 30, top.level.dir = NULL) {
  if (is.null(top.level.dir)) {
    warning("Need to supply non.null top.level.dir to AddAllScripts; no scripts generated")
    return(NULL)
  }
  AddScript(maxK = maxK, 1, file.path(top.level.dir, "sa.sa.96"))
  AddScript(maxK = maxK, 2, file.path(top.level.dir, "sp.sp"))
  AddScript(maxK = maxK, 3, file.path(top.level.dir, "sa.sa.COMPOSITE"))
  AddScript(maxK = maxK, 4, file.path(top.level.dir, "sp.sa.COMPOSITE"))
}
