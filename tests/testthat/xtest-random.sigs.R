context("Test generation of random signature profiles")

test_that("CreateRandomMutSigProfiles", {
  load("rand.syn.96.sigs.Rdata")
  set.seed(5)
  expect_equal(
    CreateRandomMutSigProfiles(
      ICAMS::catalog.row.order[["SBS96"]], 5, "prefix"),
    rand.syn.96.sigs)
})


test_that("CreateRandomSyn 5 spectra", {
  skip("fragile, diff-based test")
  skip_if_not_installed("ICAMS", minimum_version = "2.0.9")
  stopifnot(grepl("testthat", getwd(), fixed = TRUE))
  expect_true(
    # Defined in CreateRandom.R
    CreateRandomSyn(top.level.dir  = tempfile("test.random.5"),
                    seed           = 1443196,
                    num.syn.tumors = 5,
                    regress.dir    = "rdata/random.5/"))
})


test_that("CreateRandomSyn 1000 spectra", {
  skip_if_not_installed("ICAMS", minimum_version = "2.0.9")
  skip_on_cran()
  # debug(CreateOneSetOfRandomCatalogs)
  expect_true(
    CreateRandomSyn(top.level.dir  = tempfile("test.random.1000"),
                    seed           = 1443196,
                    num.syn.tumors = 1000,
                    regress.dir    = "rdata/syn.30.random.sigs/",
                    verbose        = TRUE))
})
