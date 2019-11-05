context("Test Create SBS1+SBS5 correlated datasets.")

test_that("CreateSBS1SBS5CorrelatedSyntheticData", {
  skip_if_not_installed("ICAMS", minimum_version = "2.0.9")

  tmp.top.level.dir <- "./tmp.SBS1.5/"

  expect_true(
    CreateSBS1SBS5CorrelatedSyntheticDataDemo(
      top.level.dir = "./tmp.SBS1.5/",
      regress.dir = "rdata/SBS1.5/",
      overwrite = FALSE,
      add.info = FALSE,
      unlink = TRUE)
  )
})
