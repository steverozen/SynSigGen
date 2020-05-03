test_that("ReadExposure", {
  in.exposure <- system.file("tests/test.data/tiny.exposure.csv",
                             package = "SynSigGen",
                             mustWork = TRUE)
  expect_warning(x <- ReadExposure(in.exposure))
  tfile <- tempfile()
  WriteExposure(x, tfile)
  reread.x <- ReadExposure(tfile)
  expect_equal(x, reread.x)

  x <- ReadExposure(in.exposure, check.names = FALSE)
  WriteExposure(x, tfile)
  x2 <- ReadExposure(tfile, check.names = FALSE)
  expect_equal(x, x2)

  in.exposure2 <- system.file("tests/test.data/tiny.exposure.dup.csv",
                             package = "SynSigGen",
                             mustWork = TRUE)
  expect_error(x <- ReadExposure(in.exposure2, check = FALSE))
})
