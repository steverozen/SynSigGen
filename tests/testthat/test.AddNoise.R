test_that("AddNoise", {

  in.exp <- matrix(c(1000, 2000, 2000, 2000, 4000, 1000),nrow = 3)
  rownames(in.exp) <- c("SBS1", "SBS4", "SBS13")
  colnames(in.exp) <- c("a", "b")

  reg.data <- new.env()
  load(
    system.file("tests/test.data/AddNoise.test.results.Rdata",
                package = "SynSigGen"),
    envir = reg.data)


  set.seed(20200601)
  add.noise.ret.1 <- AddNoise(input.exposure = in.exp,
                              signatures     = PCAWG7::signature$genome$SBS96)
  # Run manually if necessary to check resulting spectra:
  # PlotCatalogToPdf(ICAMS::as.catalog(retval$spectra), "foo.pdf")
  expect_equal(add.noise.ret.1, reg.data$add.noise.ret.1)

  add.noise.ret.2 <- AddNoise(input.exposure = in.exp,
                              signatures     = PCAWG7::signature$genome$SBS96,
                              n.binom.size   = 10)
  # Run manually if necessary to check resulting spectra:
  # PlotCatalogToPdf(ICAMS::as.catalog(retval2$spectra), "foo2.pdf")
  expect_equal(add.noise.ret.2, reg.data$add.noise.ret.2)
})

