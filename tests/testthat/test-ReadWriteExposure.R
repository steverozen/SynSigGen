test_that("ReadExposure", {
  in.exposure <- "testdata/tiny.exposure.csv"
  expect_warning(x <- ReadExposure(in.exposure))
  tfile <- tempfile()
  WriteExposure(x, tfile)
  reread.x <- ReadExposure(tfile)
  expect_equal(x, reread.x)

  x <- ReadExposure(in.exposure, check.names = FALSE)
  WriteExposure(x, tfile)
  x2 <- ReadExposure(tfile, check.names = FALSE)
  expect_equal(x, x2)

  in.exposure2 <- "testdata/tiny.exposure.dup.csv"
  expect_error(x <- ReadExposure(in.exposure2, check = FALSE))
})
