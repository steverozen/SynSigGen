test_that("AddNoise", {

  skip_if_not_installed("PCAWG7")

  in.exp <- matrix(c(1000, 2000, 2000, 2000, 4000, 1000),nrow = 3)
  rownames(in.exp) <- c("SBS1", "SBS4", "SBS13")
  colnames(in.exp) <- c("a", "b")

  reg.data <- new.env()
  load("testdata/AddNoise.test.results.Rdata", envir = reg.data)


  set.seed(20200601)
  add.noise.ret.1 <- AddNoise(input.exposure = in.exp,
                              signatures     = cosmicsig::COSMIC_v3.0$signature$GRCh37$SBS96,
                              n.binom.size = NULL)
  # Run manually if necessary to check resulting spectra:
  # PlotCatalogToPdf(ICAMS::as.catalog(retval$spectra), "foo.pdf")
  expect_equal(add.noise.ret.1, reg.data$add.noise.ret.1)

  add.noise.ret.2 <- AddNoise(input.exposure = in.exp,
                              signatures     = cosmicsig::COSMIC_v3.0$signature$GRCh37$SBS96,
                              n.binom.size   = 10)
  # Run manually if necessary to check resulting spectra:
  # PlotCatalogToPdf(ICAMS::as.catalog(retval2$spectra), "foo2.pdf")
  expect_equal(add.noise.ret.2, reg.data$add.noise.ret.2)
})

