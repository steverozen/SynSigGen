context("Test CreateFromReal")

test_that("CreateFromReal PancAdeno-CA", {
  skip_if_not_installed("ICAMS", minimum_version = "2.0.9")
  reg <- new.env()
  load("testdata/PancAdenoCA1000.retval.Rdata", envir = reg)

    retval <- PancAdenoCA1000(seed = 123,
                              num.syn.tumors = 5,
                              top.level.dir  = tempfile(pattern = "regress.pancadeno"),
                              regress.dir    = NULL, # "rdata/Panc-AdenoCA.123/",
                              unlink         = TRUE)
   expect_equal(retval$info.list[[1]]$sp.syn.exp,
                reg$retval$info.list[[1]]$sp.syn.exp)
})

test_that("CreateFromReal ManyTypes - Fast test - only generate 18 tumors", {
  run.manually <- FALSE
  if (!run.manually) {
    skip(
      paste("Run manully with devtools::test(filter = \"CreateFromReal\"): ",
            "excessively fragile test because of reliance on diff"))
  }
  skip_if_not_installed("ICAMS", minimum_version = "2.0.9")

  cat("\n\n=============================\n")
  cat("test_that(\"CreateFromReal Many\"\n")
  cat(RNGkind(), "\n")
  cat(.Random.seed[1:4])
  cat("\n=============================\n")

  expect_true(
  CreateFromReal(
    seed           = 191906,
    top.level.dir  = tempfile(pattern = "regress.create.many"),
    num.syn.tumors = 2,
    cancer.types   = c("Bladder-TCC",    "Eso-AdenoCA",
                       "Breast-AdenoCA", "Lung-SCC",
                       "Kidney-RCC",     "Ovary-AdenoCA",
                       "Bone-Osteosarc", "Cervix-AdenoCA",
                       "Stomach-AdenoCA"),
    sa.exp          = sa.all.real.exposures,
    sp.exp          = sp.all.real.exposures,
    regress.dir     = "rdata/many.syn.191906/",
    unlink          = TRUE))
})

test_that("CreateFromReal RCCOvary1000", {
  skip_if_not_installed("ICAMS", minimum_version = "2.0.9")
  run.manually <- FALSE
  if (!run.manually) {
    skip(
      paste("Run manully with devtools::test(filter = \"CreateFromReal\"): ",
            "excessively fragile test because of reliance on diff"))
  }
  # Last tested on Windows, 2019 08 25

  cat("\n\n=============================\n")
  cat("#1 test_that(\"CreateFromReal RCCOvary1000\"\n")
  cat(RNGkind(), "\n")
  cat(.Random.seed[1:4])
  cat("\n=============================\n")
  rkind <- RNGkind()
  RNGkind(kind = rkind[1], normal.kind = rkind[2], sample.kind = "default")
  cat("\n\n=============================\n")
  cat("#2 test_that(\"CreateFromReal RCCOvary1000\"\n")
  cat(RNGkind(), "\n")
  cat(.Random.seed[1:4])
  cat("\n=============================\n")

  tmp.top.level.dir <- "TMP.3.5.40.RCC.and.ovary"
  if (dir.exists(tmp.top.level.dir)) {
    stop("\n\n!!! Must manually delete ",
         file.path(getwd()), tmp.top.level.dir)
  }

  cat("\nAfter test, must manually delete file created by this test:\n")
  cat(file.path(getwd(), tmp.top.level.dir))
  cat("\n=============================\n")
  log <- capture_messages(
    expect_true(
      RCCOvary1000(regress.dir = "rdata/rcc.etc/",
                   top.level.dir = tmp.top.level.dir,
                   unlink = FALSE)))
  cat(log)
  # print(Sys.getenv())
})
