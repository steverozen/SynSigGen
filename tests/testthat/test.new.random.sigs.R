context("Test generation of random exposures with real and synthetic signature profiles")

test_that("Generate.SP.signatures.random.subsets()", {
  newdir <- tempfile("SP.signatures.random.subsets")
  # newdir <- "tmp.SP.signatures.random.subsets"
  # debug(GenerateMatrixRandomExpKnownSigs200)
  Generate.SP.signatures.random.subsets(top.level.dir = newdir,
                                        verbose = TRUE,
                                        overwrite = TRUE)
  log <- testthat::capture_messages(
  expect_equal(
    NewDiff4SynDataSets(
      newdir = newdir,
      # regressdirname = "rdata/SP.sig.ran.Rej/",
      regressdirname = "rdata/SP.sig.ran/",
      unlink = FALSE,      verbose = TRUE,
      long.diff = FALSE)[1]
    ,"ok"))
  cat(log)
  # print(Sys.getenv())
})
