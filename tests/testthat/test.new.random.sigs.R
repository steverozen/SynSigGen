context("Test generation of random exposures with real and synthetic signature profiles")

test_that("Generate.SP.signatures.random.subsets()", {
  newdir <- tempfile("SP.signatures.random.subsets")
  Generate.SP.signatures.random.subsets(top.level.dir = newdir)
  log <- testthat::capture_messages(
  expect_equal(
    NewDiff4SynDataSets(
      newdir = newdir,
      regressdirname = "rdata/SP.sig.ran/",
      unlink = FALSE,
      verbose = TRUE)[1]
    ,"ok"))
  cat(log)
})
