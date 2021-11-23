context("Test CreateAbstract")

test_that("Create.2.7a.7b.Abstract", {
  skip_if_not_installed("ICAMS", minimum_version = "2.0.9")
  expect_true(
    Create.2.7a.7b.Abstract(seed = 123,
                    num.syn.tumors = 5,
                    top.level.dir  = tempfile(pattern = "regress.2.7a.7b.abstract"),
                    regress.dir    = "rdata/abst.2.7a.7b/",
                    unlink         = TRUE))
})

