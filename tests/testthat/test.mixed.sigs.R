context("Test CreateFromReal")

test_that("CreateFromReal PancAdeno-CA", {
  skip_if_not_installed("ICAMS", minimum_version = "2.0.9")
  expect_true(
    PancAdenoCA1000(seed = 123,
                    num.syn.tumors = 5,
                    top.level.dir  = tempfile(pattern = "regress.pancadeno"),
                    regress.dir    = "rdata/Panc-AdenoCA.123/",
                    unlink         = TRUE))
})

test_that("CreateFromReal Many", {
  skip_if_not_installed("ICAMS", minimum_version = "2.0.9")
  if (FALSE) {
  cat("\n\n=============================\n")
  cat("test_that(\"CreateFromReal Many\"\n")
  cat(RNGkind(), collapse = " ")
  cat("\n")
  cat(.Random.seed[1:4], colapse = " ")
  cat("\n=============================\n")
  }
  expect_true(
  CreateFromReal(
    seed           = 123,
    top.level.dir  = tempfile(pattern = "regress.create.many"),
    num.syn.tumors = 2,
    cancer.types   = c("Bladder-TCC",    "Eso-AdenoCA",
                       "Breast-AdenoCA", "Lung-SCC",
                       "Kidney-RCC",     "Ovary-AdenoCA",
                       "Bone-Osteosarc", "Cervix-AdenoCA",
                       "Stomach-AdenoCA"),
    sa.exp          = sa.all.real.exposures,
    sp.exp          = sp.all.real.exposures,
    regress.dir     = "rdata/many.syn.123/",
    unlink          = TRUE))
})

test_that("CreateFromReal RCCOvary1000", {
  skip_if_not_installed("ICAMS", minimum_version = "2.0.9")
  # skip("Takes too long")
  if (FALSE) {
  cat("\n\n=============================\n")
  cat("#1 test_that(\"CreateFromReal RCCOvary1000\"\n")
  cat(RNGkind(), collapse = " ")
  cat("\n")
  cat(.Random.seed[1:4], colapse = " ")
  cat("\n=============================\n")
  rkind <- RNGkind()
  RNGkind(kind = rkind[1], normal.kind = rkind[2], sample.kind = "default")
  cat("\n\n=============================\n")
  cat("#2 test_that(\"CreateFromReal RCCOvary1000\"\n")
  cat(RNGkind(), collapse = " ")
  cat("\n")
  cat(.Random.seed[1:4], colapse = " ")
  cat("\n=============================\n")
  }

  log <- capture_messages(

    expect_true(
      RCCOvary1000(regress.dir = "rdata/rcc.etc/", unlink = FALSE)))
  cat(log)
  # print(Sys.getenv())
})
