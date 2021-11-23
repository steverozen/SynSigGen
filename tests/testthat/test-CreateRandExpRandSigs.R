context("Generate spectra with random exposures and random synthetic signature profiles")

test_that("GenerateAllRandomSP", {
  run.manually <- FALSE # If FALSE then skipped
  if (run.manually) {
    newdir <- "tmp.SP.rand.exp.rand.sigs"
    cat("\nremember to manually delete", newdir, "\n")
    unlink <- FALSE
  } else {
    skip(
      paste("Run manully with devtools::test(filter = \"RandExpRandSig\"): ",
            "excessively fragile test because of reliance on diff"))
  }

  GenerateAllRandomSP(
    top.level.dir  = newdir,
    verbose        = TRUE,
    overwrite      = TRUE,
    unlink         = unlink,
    num.replicates = 1)
  unlink.ret <- unlink(file.path(newdir, "log.txt")) # Remove from the diff
  if (unlink.ret != 0) cat("failed to unlink log file")
  log <- testthat::capture_messages(
  expect_equal(
    NewDiff4SynDataSets(
      newdir = newdir,
      regressdirname = "rdata/r.exp.r.sigs/",
      unlink = FALSE,      verbose = TRUE,
      long.diff = FALSE)[1]
    ,"ok"))
  cat(log)
})
