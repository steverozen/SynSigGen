#' Create a test data set based on >= 1 tumor types.
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
#' @param sa.exp SignatureAnalyzer exposures from which to select cancer types
#'        specified by \code{cancer.type.strings}. In the
#'        column names of \code{sa.exp} the cancer type string
#'        should be separated from the sample identifier by two colons
#'        (::).
#'
#' @param sp.exp SigProfiler exposures from which to select cancer types
#'        specified by \code{cancer.type.strings}. In the
#'        column names of \code{sp.exp} the cancer type string
#'        should be separated from the sample identifier by two colons
#'        (::).
#
#' @param verbose If > 0, cat various messages.
#'
#' @param bladder.regress.hack For use by \code{\link{BladderSkin1000}}.
#'     Forces use of non-hyper-mutated exposures for bladder-TCC even if
#'     \code{sa.exp} and \code{sp.exp} include hyper-mutated exposures.
#'
#' @export

CreateMixedTumorTypeSyntheticData <- function(top.level.dir,
                                              cancer.type.strings,
                                              num.syn.tumors,
                                              overwrite = FALSE,
                                              sa.exp = sa.all.real.exposures,
                                              sp.exp = sp.all.real.exposures,
                                              verbose = FALSE,
                                              bladder.regress.hack = FALSE) {

    odigits <- getOption("digits")
    options(digits = 9)
    on.exit(options(digits = odigits))

    info.list <-
      lapply(cancer.type.strings,
             function(ca.type.str) {
               if (verbose) message("\n\nProcessing", ca.type.str, "\n\n\n")
               local.sa.exp <- sa.exp
               local.sp.exp <- sp.exp
               if (bladder.regress.hack && ca.type.str == "Bladder-TCC") {
                 local.sa.exp <- SynSigGen::sa.no.hyper.real.exposures
                 local.sp.exp <- SynSigGen::sp.no.hyper.real.exposures
                 message("bladder.regress.hack deployed for ", ca.type.str)
               }
               retval <-
                 SAAndSPSynDataOneCAType(
                   sa.real.exp    = local.sa.exp,
                   sp.real.exp    = local.sp.exp,
                   ca.type        = ca.type.str,
                   num.syn.tumors = num.syn.tumors,
                   file.prefix    = ca.type.str,
                   top.level.dir  = top.level.dir)
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
      dir = NULL, # "sa.sa.COMPOSITE",
      WriteCatCOMPOSITE,
      overwrite = overwrite,
      my.dir = file.path(top.level.dir, "sa.sa.COMPOSITE"))

    CreateAndWriteCatalog(
      sa.96.sigs,
      sa.exp,
      dir = NULL, # "sa.sa.96",
      WriteCatalog,
      overwrite = overwrite,
      my.dir = file.path(top.level.dir, "sa.sa.96"))

    # Create catalogs of synthetic mutational spectra
    # based on SigProfiler attributions

    # First we need the matching between SigProfiler and
    # SignatureAnalyzers signatures.

    sp.sa.map.info <-
      MapSPToSASignatureNamesInExposure(sp.exp)

    CreateAndWriteCatalog(
      sa.COMPOSITE.sigs,
      sp.sa.map.info$exp2,
      dir = NULL, # "sp.sa.COMPOSITE",
      WriteCatCOMPOSITE,
      overwrite = overwrite,
      my.dir = file.path(top.level.dir, "sp.sa.COMPOSITE"))


    if (verbose) print(sp.sa.map.info$sp.to.sa.sig.match)

    CreateAndWriteCatalog(
      sp.sigs,
      sp.exp,
      dir = NULL, # "sp.sp",
      WriteCatalog,
      overwrite = overwrite,
      my.dir = file.path(top.level.dir, "sp.sp"))

    # AddAllScripts(maxK = 50)

    invisible(list(info.list = info.list,
                   sp.sa.map.info =  sp.sa.map.info))

  }


RNGMessages <- function(prefix = NULL, message.fn = message) {
  message.fn("\n")
  if (!is.null(prefix)) message.fn(prefix, "\n")
  message.fn("RNGkind = ", paste(RNGkind(), collapse = " "), "\n")
  message.fn(".Random.seed[1:4] = ",
             paste(.Random.seed[1:4], colapse = " "),
             "\n")
}


#' Create a specific synthetic data set based on real exposures in one or more cancer types.
#'
#' Create a full SignatureAnalyzer / SigProfiler test data set for a
#' set of various tumor types.
#'
#' @param seed A random seed to use.
#'
#' @param top.level.dir The directory in which to put the output; will
#'        be created if necessary.
#'
#' @param enclosing.dir Deprecated; create the output in a subdirectory of this directory.
#'
#' @param num.syn.tumors The number of tumors to create \strong{for each cancer
#'    type} in \code{cancer.types}.
#'
#' @param cancer.types Search \code{sa.exp} and \code{sp.exp}
#' for exposures from tumors matching
#' these strings. Each string should identify one tumor type, for
#' some definition of tumor type. Probably the tumors in each type
#' should be non-overlapping, but the code does not enforce this and
#' does not care.
#'
#' @param data.suite.name Deprecated; the directory created will be
#'   \code{file.path(enclosing.dir, paste0(data.suite.name, ".", seed))}.
#'
#' @param sa.exp A matrix of exposures; this function will use the
#'        columns with column names beginning \code{paste0(cancer.type, "::")}.
#'
#' @param sp.exp A matrix of exposures; this function will use the
#'        columns with column names beginning \code{paste0(cancer.type, "::")}.
#'
#' @param overwrite If TRUE, overwrite existing directories and files.
#'
#' @param regress.dir If not \code{NULL}, compare the result to
#'    the contents of this directory with a \code{diff}.
#'
#' @param unlink If \code{TRUE} and \code{!is.null(regress.dir)}, then
#'       unlink the result directory if there are no differences.
#'
#' @param verbose If \code{TRUE} print various informative messages.
#'
#' @param bladder.regress.hack Set this to \code{TRUE} to handle
#'        mixed "all" and "no hyper" signature sets for the
#'        regression test for \code{\link{BladderSkin1000}}.
#'
#' @export

CreateFromReal <- function(seed,
                           top.level.dir   = NULL,
                           enclosing.dir   = NULL,
                           num.syn.tumors,
                           cancer.types,
                           data.suite.name = NULL,
                           sa.exp          = SynSigGen::sa.all.real.exposures,
                           sp.exp          = SynSigGen::sp.all.real.exposures,
                           overwrite       = TRUE,
                           regress.dir     = NULL,
                           unlink          = FALSE,
                           verbose         = FALSE,
                           bladder.regress.hack = FALSE) {



  if (verbose) RNGMessages("In CreateFromReal before set.seed", cat)
    # rkind <- RNGkind()
  # RNGkind(kind = rkind[1], normal.kind = rkind[2], sample.kind = "default")
  suppressWarnings(set.seed(seed, sample.kind = "Rounding"))
  if (verbose) RNGMessages("In CreateFromReal after set.seed", cat)

  if (!is.null(top.level.dir)) {
    if (!is.null(enclosing.dir))   stop("Do not specify top.level.dir and enclosing.dir")
    if (!is.null(data.suite.name)) stop("Do not specify data.suite.name and enclosing.dir")
  } else {
    if (is.null(enclosing.dir)) {
      stop("If top.level.dir is NULL enclosing.dir must be non-NULL")
    }
    if (is.null(data.suite.name)) {
      stop("If top.level.dir is NULL data.suite.name must be non-NULL")
    }
    top.level.dir <-
      file.path(enclosing.dir, paste0(data.suite.name, ".", seed))
  }

  MustCreateDir(top.level.dir, overwrite)

  retval <-
    CreateMixedTumorTypeSyntheticData(
      top.level.dir = top.level.dir,
      cancer.type.strings = cancer.types,
      num.syn.tumors = num.syn.tumors,
      sa.exp = sa.exp,
      sp.exp = sp.exp,
      overwrite = overwrite,
      verbose   = verbose,
      bladder.regress.hack = bladder.regress.hack
    )

  if (!is.null(regress.dir)) {
    diff.result <-
      NewDiff4SynDataSets(newdir         = top.level.dir,
                          regressdirname = regress.dir,
                          unlink         = unlink,
                          verbose        = verbose,
                          long.diff      = FALSE)
    return(diff.result[1] == "ok")
  }
  return(retval)
}

#' Create 1000 synthetic pancreatic adenocarcinoma spectra.
#'
#' This function generates synthetic tumor spectra with mutational
#' signature prevalence and mutation load similar to pancreatic
#' adenocarcinoma in PCAWG cohort.
#'
#' This function replaces \code{data-raw/Create.pancreas.Rmd} in GitHub
#' repository \code{steverozen/SynSig}. With default arguments, this
#' function generates the same results as \code{data-raw/Create.pancreas.Rmd}.
#'
#' Data set generated by this function can be found at Synapse with Synapse ID:
#' \href{https://www.synapse.org/#!Synapse:syn18500212}{syn18500212}.
#'
#' @param seed A random seed to use.
#'
#' @param regress.dir If not \code{NULL}, compare the result to
#' the contents of this directory with a \code{diff}.
#'
#' @param num.syn.tumors Generate this number of synthetic tumors.
#'
#' @param top.level.dir The directory in which to put the output; will
#' be created if necessary.
#'
#' @param unlink If \code{TRUE} and \code{!is.null(regress.dir)}, then
#' unlink the result directory if there are no differences.
#'
#' @export

PancAdenoCA1000 <- function(
  seed           = 191907,
  regress.dir    = "data-raw/long.test.regression.data/syn.pancreas/",
  num.syn.tumors = 1000,
  top.level.dir  = "../Pan-AdenoCA",
  unlink         = FALSE) {
  CreateFromReal(
    seed           = seed,
    top.level.dir  = top.level.dir,
    num.syn.tumors = num.syn.tumors,
    cancer.types   = "Panc-AdenoCA",
    sa.exp         = SynSigGen::sa.no.hyper.real.exposures,
    sp.exp         = SynSigGen::sp.no.hyper.real.exposures,
    regress.dir    = regress.dir,
    unlink         = unlink
  )
}


#' Create synthetic spectra based on renal cell carcinoma and ovarian adenocarcinoma
#'
#' Creates spectra dataset consists of 500 synthetic renal cell carcinoma (RCC)
#' with high prevalence and mutation load from SBS5 and SBS40 signatures,
#' and 500 synthetic ovarian adenocarcinoma with high prevalence and
#' mutation load from SBS3. This dataset challenges the computational approaches
#' as these three signatures are "flat" signatures hard to be extracted accurately.
#'
#' This function Replaces the first part of \code{data-raw/Create.3.5.40.Rmd} in
#' GitHub repository \code{steverozen/SynSig}. With default arguments, this
#' function generates the same results as the first part of
#' \code{data-raw/Create.3.5.40.Rmd}.
#'
#' The second half of \code{data-raw/Create.3.5.40.Rmd} in \code{steverozen/SynSig}
#' is replaced by \code{\link{Create.3.5.40.Abstract}}.
#'
#' Data set generated by this function can be found at Synapse with Synapse ID:
#' \href{https://www.synapse.org/#!Synapse:syn18500214}{syn18500214}.
#'
#' @param seed A random seed to use.
#'
#' @param unlink The directory in which to put the output; will
#' be created if necessary.
#'
#' @param regress.dir If not \code{NULL}, compare the result to
#'    the contents of this directory with a \code{diff}.
#'
#' @param top.level.dir The directory in which to put the output; will
#' be created if necessary.
#'
#'
#' @export
RCCOvary1000 <- function(seed = 191905,
                         unlink = FALSE,
                         regress.dir = NULL,
                         top.level.dir = "tmp.3.5.40.RCC.and.ovary") {
  CreateFromReal(
    seed           = seed,
    num.syn.tumors = 500,
    cancer.types   = c("Kidney-RCC", "Ovary-AdenoCA" ),
    top.level.dir  = top.level.dir,
    sa.exp         = SynSigGen::sa.no.hyper.real.exposures,
    sp.exp         = SynSigGen::sp.no.hyper.real.exposures,
    regress.dir    = regress.dir,
    unlink         = unlink)
}


#' Generate synthetic data sets modeled on bladder TCC and skin melanoma.
#'
#' Creates spectra dataset consists of 500 synthetic bladder
#' transitional cell carcinoma with high prevalence and mutation
#' load from SBS2, and 500 synthetic skin melanoma
#' with high prevalence and mutation load from SBS7a and SBS7b. This
#' dataset challenges the computational approaches as SBS2 has a similar
#' pattern to the mixture of SBS7a and SBS7b, thus the existence of these
#' signatures may interfere computational approaches from accurately
#' extracting these signatures.
#'
#' This function replaces the first part of \code{data-raw/Create.2.7a.7b.Rmd}
#' in GitHub repository \code{steverozen/SynSig}. With default arguments, this
#' function generates the same results as the first part of
#' \code{data-raw/Create.2.7a.7b.Rmd}.
#'
#' #' The second half of \code{data-raw/Create.2.7a.7b.Rmd}
#' is replaced by \code{\link{Create.2.7a.7b.Abstract}}.
#'
#' Data set generated by this function can be found at Synapse with Synapse ID:
#' \href{https://www.synapse.org/#!Synapse:syn18500217}{syn18500217}.
#'
#' @param seed A random seed to use.
#'
#' @param regress Whether to compare the result with
#' local copy of dataset using a \code{diff}.
#'
#' @export

BladderSkin1000 <- function(seed = 191906, regress = FALSE) {
  if (regress) {
    regress.dir <-
      "data-raw/long.test.regression.data/syn.2.7a.7b.bladder.and.melanoma/"
  } else regress.dir <- NULL

  CreateFromReal(
    seed            = seed,
    enclosing.dir   = "..",
    num.syn.tumors  = 500,
    cancer.types    = c("Bladder-TCC", "Skin-Melanoma" ),
    data.suite.name = "2.7a.7b.bladder.and.melanoma",
    sa.exp          = SynSigGen::sa.all.real.exposures,
    sp.exp          = SynSigGen::sp.all.real.exposures,
    regress.dir = regress.dir,
    bladder.regress.hack = TRUE
  )
}

#' Create a specific synthetic data set of 2,700 tumors.
#'
#'
#' Data set generated by this function can be found at Synapse with Synapse ID:
#' \href{https://www.synapse.org/#!Synapse:syn18500213}{syn18500213}.
#'
#' @param seed A random seed to use.
#'
#' @param regress Whether to compare the result with
#' local copy of dataset using a \code{diff}.
#'
#' @keywords export

ManyTypes2700 <- function(seed = 191906, regress = FALSE) {
  if (regress) {
    regress.dir <-
      "data-raw/long.test.regression.data/syn.many.types/"
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
    PancAdenoCA1000(seed = seed, regress.dir = NULL)
  }
}

three.5.40.for.ludmil.2019.08.17 <- function() {
  seeds <- sample(10000, size = 9)
  for (seed in seeds) {
    RCCOvary1000(seed = seed)
  }
}

bladder.skin.for.ludmil.2019.08.18 <- function() {
  seeds <- sample(10000, size = 9)
  for (seed in seeds) {
    BladderSkin1000(seed = seed)
  }
}
