#' diff new directory / files against regression data for testing.
#'
#' @param dirname the root name of the directories to diff.
#'
#' @param unlink if \code{TRUE} unlink \code{tmpdirname}, but do not unlink
#' if there are diffs.
#'
#' @return The output of the diff command.
#'
#' @export

Diff4SynDataSets <- function(dirname, unlink) {
  regressdirname <- paste0("data-raw/long.test.regression.data/", dirname)
  tmpdirname <- paste0("tmp.", dirname)
  return(NewDiff4SynDataSets(newdir = tmpdirname,
                             regressdirname = regressdirname,
                             unlink = unlink))
}

#' Diff two directories or files.
#'
#' @param newdir the path of \code{dir2} for a folder to be recursively compared with
#' \code{dir1}; it can also be the path of a single file \code{file2} to diff with
#' \code{file1}.
#'
#' @param regressdirname the path of \code{dir2} for a folder to be recursively compared with
#' \code{dir1}; it can also be the path of a single file \code{file2} to diff with
#' \code{file1}.
#'
#' @param unlink if \code{TRUE} unlink \code{newdir}, but do not unlink
#' if there are diffs.
#'
#' @param verbose Whether to display additional R messages.
#'
#' @param long.diff If \code{TRUE}, invoke \code{diff -r} (detailed text information
#' even if the two files/folders are the same); if \code{FALSE}, invoke \code{diff -rq}
#' (detailed text information only if two files/folders are different).
#' (Default: FALSE)
#'
#' @return The output of the diff command.
#'
#' @export
NewDiff4SynDataSets <-
  function(newdir, regressdirname, unlink, verbose = FALSE, long.diff = FALSE) {

  if (!dir.exists(regressdirname)) {
    stop("regressdirname ", regressdirname, " does not exist\n",
         "getwd() = ", getwd())
  }
  if (!dir.exists(newdir)) {
    stop("new directory ", newdir, " does not exist\n",
         "getwd() = ", getwd())
  }
  if (long.diff) {
    flags <- "-r"
  } else {
    flags <- "-rq"
  }
  cmd.result <-
    system2("diff", c(flags, newdir, regressdirname),
            stderr = TRUE, stdout = TRUE) # Capture all output
  if (length(cmd.result) == 0) {
    # No differences
    if (unlink) {
      unlink.res <- unlink(newdir, recursive = TRUE, force = TRUE)
      if (unlink.res != 0) {
        warning("failed to unlink ", newdir)
        return("failed to unlink")
      }
    }
    if (verbose) message("\nok\n")
    return("ok")
  }

  cmd.result <-
    c("diff", paste("diff", flags, newdir, regressdirname), cmd.result)
  message("\nThere was a difference, investigate\n",
          paste0(cmd.result, "\n"))
  return(cmd.result)
}
