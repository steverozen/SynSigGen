context("Generate spectra with random exposures with real signature profiles")

test_that("Generate.SP.signatures.random.subsets()", {
  run.manually <- FALSE # If FALSE then skipped
  if (run.manually) {
    newdir <- "tmp.SP.signatures.random.subsets"
    cat("remember to manually delete", newdir, "\n")
    unlink <- FALSE
  } else {
    skip(
      paste("Run manully with devtools::test(filter = \"RandExpKnownSig\"): ",
            "excessively fragile test because of reliance on diff"))
  }

  # debug(GenerateMatrixRandomExpKnownSigs200)
  Generate.SP.signatures.random.subsets(top.level.dir  = newdir,
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
      regressdirname = "rdata/r.exp.k.sigs/",
      unlink = FALSE,      verbose = TRUE,
      long.diff = FALSE)[1]
    ,"ok"))
  cat(log)
})
