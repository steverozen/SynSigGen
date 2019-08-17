#' Create a full SignatureAnalyzer / SigProfiler test data set for a
#' set of various tumor types.
#'
#' @param top.level.dir Path to top level of directory structure to be created.
#'
#' @param cancer.type.strings Search the PCAWG data for tumors matching
#' these strings. Each string should identify one tumor type, for
#' some definition of tumor type. Probably the tumors in each type
#' should be non-overlapping, but the code does not enforce this and
#' does not care.
#'
#' @param num.syn.tumors Number of synthetic tumors to create
#' for each cancer type.
#'
#' @param overwrite If TRUE, overwrite existing directories / files.
#'
#' @param verbose If > 0, cat various messages.
#'
#' @export

CreateMixedTumorTypeSyntheticData <- function(top.level.dir,
                                              cancer.type.strings,
                                              num.syn.tumors,
                                              overwrite = FALSE,
                                              sa.exp = sa.all.real.exposures,
                                              sp.exp = sp.all.real.exposures,
                                              verbose = 0) {

    SetNewOutDir(top.level.dir, overwrite)

    odigits <- getOption("digits")
    options(digits = 9)
    on.exit(options(digits = odigits))

    info.list <-
      lapply(cancer.type.strings,
             function(ca.type.str) {
               if (verbose) cat("\n\nProcessing", ca.type.str, "\n\n\n")
               retval <-
                 SAAndSPSynDataOneCAType(
                   sa.real.exp = sa.exp,
                   sp.real.exp = sp.exp,
                   ca.type = ca.type.str,
                   num.syn.tumors = num.syn.tumors,
                   file.prefix = ca.type.str,
                   top.level.dir = top.level.dir)
               return(retval)
             })

    # info.list has both SignatureAnalyzer and
    # SigProfiler exposures.

    sa.exposures <- lapply(info.list, function(x) x$sa.syn.exp)
    sa.exp <- MergeExposures(sa.exposures)
    sa.exp <- sa.exp[rowSums(sa.exp) > 0.5, ]
    if (verbose) message("Dimension sa.exp", dim(sa.exp))

    sp.exposures <- lapply(info.list, function(x) x$sp.syn.exp)
    sp.exp <- MergeExposures(sp.exposures)
    sp.exp <- sp.exp[rowSums(sp.exp) > 0.5, ]
    if (verbose) message("Dimension sp.exp", dim(sp.exp))


    # We will need the exposures later when evaluating the attributed signatures
    WriteExposure(sa.exp, file.path(top.level.dir, "sa.exposure.csv"))
    WriteExposure(sp.exp, file.path(top.level.dir, "sp.exposure.csv"))

    # Create catalogs of synthetic mutational spectra
    # based on SignatureAnalyzer attributions

    CreateAndWriteCatalog(
      sa.COMPOSITE.sigs,
      sa.exp,
      "sa.sa.COMPOSITE",
      WriteCatCOMPOSITE,
      overwrite = overwrite)

    CreateAndWriteCatalog(
      sa.96.sigs,
      sa.exp,
      "sa.sa.96",
      WriteCatalog,
      overwrite = overwrite)

    # Create catalogs of synthetic mutational spectra
    # based on SigProfiler attributions

    # First we need the matching between SigProfiler and
    # SignatureAnalyzers signatures.

    sp.sa.map.info <-
      MapSPToSASignatureNamesInExposure(sp.exp)

    CreateAndWriteCatalog(
      sa.COMPOSITE.sigs,
      sp.sa.map.info$exp2,
      "sp.sa.COMPOSITE",
      WriteCatCOMPOSITE,
      overwrite = overwrite)


    if (verbose) print(sp.sa.map.info$sp.to.sa.sig.match)

    CreateAndWriteCatalog(
      sp.sigs,
      sp.exp,
      "sp.sp",
      WriteCatalog,
      overwrite = overwrite)

    # AddAllScripts(maxK = 50)

    invisible(list(info.list = info.list,
                   sp.sa.map.info =  sp.sa.map.info))

  }

#' Create a specific synthetic data set of 2,700 tumors.
#'
#' @param regress If \code{TRUE}, then compare to data in \code{data-raw}
#' and report any differences; if no differences, unlink the result
#' directory.
#'
#' @export

Create.syn.many.types <- function(regress = FALSE, seed = NULL, unlink = FALSE) {
  if (is.null(seed)) {
    suppressWarnings(RNGkind(sample.kind = "Rounding"))
    # For compatibility with R < 3.6.0
    set.seed(191906)
    top.level.dir <- "tmp.syn.many.types"

  } else {
    set.seed(seed)
    top.level.dir <- paste0("../2700.tumors.seed.", seed)
  }
  num.syn.tumors <- 300 # number of tumor of each type

  cancer.types <- c("Bladder-TCC",    "Eso-AdenoCA",
                    "Breast-AdenoCA", "Lung-SCC",
                    "Kidney-RCC",     "Ovary-AdenoCA",
                    "Bone-Osteosarc", "Cervix-AdenoCA",
                    "Stomach-AdenoCA")
  retval <-
    CreateMixedTumorTypeSyntheticData(
      top.level.dir = top.level.dir,
      cancer.type.strings = cancer.types,
      num.syn.tumors = num.syn.tumors,
      overwrite = TRUE
    )

  if (regress) {
    diff.result <- Diff4SynDataSets("syn.many.types", unlink = unlink)
    if (diff.result[1] != "ok") {
      message("\nThere was a difference, investigate\n",
              paste0(diff.result, "\n"))
    } else {
      message("\nok\n")
    }
  }

  invisible(retval)
}

for.ludmil.2019.08.16 <- function() {
  seeds <- sample(10000, size = 9)
  for (seed in seeds) {
    Create.syn.many.types(seed = seed)
  }
}

#' Create a specific synthetic data set of 1000 pancreatic adenocarcinomas.
#'
#' @param regress If \code{TRUE}, then compare to data in \code{data-raw}
#' and report any differences; if no differences, unlink the result
#' directory.
#'
#' @export


#' Create a specific synthetic data set based on real exposures in one or more cancer types.
#'
#' @keywords internal

CreateFromReal <- function(seed,
                           enclosing.dir,
                           num.syn.tumors,
                           cancer.types,
                           data.suite.name,
                           sa.exp      = sa.all.real.exposures,
                           sp.exp      = sp.all.real.exposures,
                           overwrite   = TRUE,
                           regress.dir = NULL,
                           unlink      = FALSE) {
  set.seed(seed)

  top.level.dir <-
    file.path(enclosing.dir, paste0(data.suite.name, ".", seed))

  retval <-
    CreateMixedTumorTypeSyntheticData(
      top.level.dir = top.level.dir,
      cancer.type.strings = cancer.types,
      num.syn.tumors = num.syn.tumors,
      sa.exp = sa.exp,
      sp.exp = sp.exp,
      overwrite = overwrite
    )

  if (!is.null(regress.dir)) {
    diff.result <-
      NewDiff4SynDataSets(top.level.dir, regress.dir, unlink = unlink)
    if (diff.result[1] != "ok") {
      message("\nThere was a difference, investigate\n",
              paste0(diff.result, "\n"))
    } else {
      message("\nok\n")
    }
  }

  invisible(retval)
}


PancAdenoCA1000 <- function(seed = 191907, regress = FALSE) {
  if (regress) {
    regress.dir <- "data-raw/long.test.regression.data/syn.pancreas/"
  } else regress.dir <- NULL
  CreateFromReal(
    seed           = seed,
    enclosing.dir = "..",
    num.syn.tumors = 1000,
    cancer.types   = "Panc-AdenoCA",
    data.suite.name = "Panc-AdenoCA",
    sa.exp      = sa.no.hyper.real.exposures,
    sp.exp      = sp.no.hyper.real.exposures,
    regress.dir = regress.dir
  )
}

RCCOvary1000 <- function(seed = 191905, regress = FALSE) {
  if (regress) {
    regress.dir <- "data-raw/long.test.regression.data/syn.3.5.40.rcc.and.ovary/"
  } else regress.dir <- NULL
  CreateFromReal(
    seed           = seed,
    enclosing.dir = "..",
    num.syn.tumors = 500,
    cancer.types   = c("Kidney-RCC", "Ovary-AdenoCA" ),
    data.suite.name = "3.5.40.RCC.and.ovary",
    sa.exp      = sa.no.hyper.real.exposures,
    sp.exp      = sp.no.hyper.real.exposures,
    regress.dir = regress.dir
  )
}


ManyTypes2700 <- function(seed = 191906, regress = FALSE) {
  if (regress) {
    regress.dir <- "data-raw/long.test.regression.data/syn.many.types/"
  } else regress.dir <- NULL
  # suppressWarnings(RNGkind(sample.kind = "Rounding"))
  # For compatibility with R < 3.6.0
  CreateFromReal(
    seed           = seed,
    enclosing.dir  = "..",
    num.syn.tumors = 300,
    cancer.types   = c("Bladder-TCC",    "Eso-AdenoCA",
                       "Breast-AdenoCA", "Lung-SCC",
                       "Kidney-RCC",     "Ovary-AdenoCA",
                       "Bone-Osteosarc", "Cervix-AdenoCA",
                       "Stomach-AdenoCA"),
    data.suite.name = "Many.types",
    sa.exp          = sa.all.real.exposures,
    sp.exp          = sp.all.real.exposures,
    regress.dir     = regress.dir
  )

}

for.ludmil.2019.08.17 <- function() {
  seeds <- sample(10000, size = 9)
  for (seed in seeds) {
    ManyTypes2700(seed = seed)
  }
}

panc.for.ludmil.2019.08.17 <- function() {
  seeds <- sample(10000, size = 9)
  for (seed in seeds) {
    PancAdenoCA1000(seed = seed)
  }
}

three.5.40.for.ludmil.2019.08.17 <- function() {
  seeds <- sample(10000, size = 9)
  for (seed in seeds) {
    RCCOvary1000(seed = seed)
  }
}

