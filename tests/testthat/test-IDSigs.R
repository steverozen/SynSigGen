test_that("Synthetic exposures for indels", {
  skip_if_not_installed("ICAMS", minimum_version = "2.0.9")


TestIDSigParams <- function() {
  retval <- data.frame(
    ID1 = c(prob = 2590 / 2780, mean = 1.9, stdev = 0.7),
    ID2 = c(prob = 1935 / 2780, mean = 1.7, stdev = 0.8),
    ID5 = c(prob = 1538 / 2780, mean = 2.0, stdev = 0.5))
  return(retval)
}

num.samples <- 200

param.partial <- TestIDSigParams()
set.seed(20200601)
expos.partial <- GenerateSyntheticExposures(param.partial,
                                            num.samples = num.samples,
                                            name = "Syn.indel")

param.prob1 <- param.partial
param.prob1["prob", ] <- 1 # More challenging data set

set.seed(20200601)
expos.prob1 <- GenerateSyntheticExposures(param.prob1,
                                          num.samples = num.samples,
                                          name = "Syn.indel")
partial.dir <- tempfile()
prob1.dir   <- tempfile()
dir.create(partial.dir)
dir.create(prob1.dir)

write.csv(param.partial, file = file.path(partial.dir, "param.partial.csv"))
write.csv(param.prob1,   file = file.path(prob1.dir,   "param.prob1.csv"))

ret.partial <- CreateAndWriteCatalog(
  sigs         = PCAWG7::signature$genome$ID,
  exp          = expos.partial,
  dir          = NULL,
  write.cat.fn = ICAMS::WriteCatalog,
  overwrite    = TRUE,
  my.dir       = partial.dir)

ret.prob1 <- CreateAndWriteCatalog(
  sigs         = PCAWG7::signature$genome$ID,
  exp          = expos.prob1,
  dir          = NULL,
  write.cat.fn = ICAMS::WriteCatalog,
  overwrite    = TRUE,
  my.dir       = prob1.dir)

# This is how the regression data were generated:
#
# save(expos.partial, expos.prob1, ret.partial, ret.prob1,
#     file = "tests/test.data/ID.test.results.Rdata")

reg.data <- new.env()
load(
  system.file("tests/test.data/ID.test.results.Rdata",
              package = "SynSigGen"),
  envir = reg.data)

testthat::expect_equal(expos.prob1,   reg.data$expos.prob1)
testthat::expect_equal(expos.partial, reg.data$expos.partial)
testthat::expect_equal(ret.prob1,     reg.data$ret.prob1)
testthat::expect_equal(ret.partial,   reg.data$ret.partial)
})
