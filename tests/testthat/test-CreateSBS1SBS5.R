context("Test Create SBS1+SBS5 correlated datasets.")

test_that("CreateSBS1SBS5CorrelatedSyntheticData", {

  run.manually <- FALSE # If FALSE then skipped
  if (run.manually) {
    skip_if_not_installed("ICAMS", minimum_version = "2.0.9")

    tmp.top.level.dir <- "./tmp.SBS1.5/"

    expect_true(
      CreateSBS1SBS5CorrelatedSyntheticData(
        top.level.dir = "./tmp.SBS1.5/",
        regress.dir = "rdata/SBS1.5/",
        overwrite = FALSE,
        add.info = FALSE,
        unlink = TRUE)
    )
  } else {
    skip(
      paste0("Run manully with devtools::test(filter = \"CreateSBS1SBS5\"): ",
            "excessively fragile test because of reliance on diff. \n",
            "Generated exposures of different system may differ at the 12th decimal place, ",
            "yet it does not affect the values for R core because this difference is smaller than eps."))

  }
})
