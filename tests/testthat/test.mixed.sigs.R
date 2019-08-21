context("Test CreateFromReal")

test_that("CreateFromReal PancAdeno-CA", {
  expect_true(
    PancAdenoCA1000(seed = 123,
                    num.syn.tumors = 5,
                    regress.dir = "rdata/Panc-AdenoCA.123/"))
})

test_that("CreateFromReal Many", {
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
    unlink          = TRUE)
})

test_that("CreateFromReal RCCOvary1000", {
  expect_identical(
    RCCOvary1000(regress.dir = "rdata/rcc.etc/", unlink = TRUE),
    TRUE)
})
