context("Test CreateFromReal")

test_that("CreateFromReal RCCOvary1000", {
  expect_identical(
    RCCOvary1000(regress.dir = "rdata/rcc.etc/", unlink = TRUE),
    TRUE)
})
