three.4.40.abstract.for.ludmil.2019.08.18 <- function() {
  seeds <- sample(10000, size = 9)
  for (seed in seeds) {
    Create.3.4.40.Abstract(seed = seed, overwrite = FALSE, regress.dir = NULL)
  }
}


Create.3.4.40.Abstract <- function(seed        = 44,
                                   overwrite   = TRUE,
                                   regress.dir = "data-raw/long.test.regression.data/syn.3.5.40.abst/") {

  top.level.dir <- paste0("../syn.3.5.40.abst.", seed)

  set.seed(seed)
  num.syn.tumors <- 1000

  sa.kidney.exp <- GetExpForOneCancerType("Kidney-RCC",
                                          SynSigGen::sa.no.hyper.real.exposures)
  sa.kidney.parms <- GetSynSigParamsFromExposures(sa.kidney.exp)
  sa.ovary.exp <- GetExpForOneCancerType("Ovary-AdenoCA",
                                         SynSigGen::sa.no.hyper.real.exposures)
  sa.ovary.parms  <- GetSynSigParamsFromExposures(sa.ovary.exp)


  sp.kidney.exp <- GetExpForOneCancerType("Kidney-RCC",
                                          SynSigGen::sp.no.hyper.real.exposures)
  sp.kidney.parms <- GetSynSigParamsFromExposures(sp.kidney.exp)
  sp.ovary.exp <- GetExpForOneCancerType("Ovary-AdenoCA",
                                         SynSigGen::sp.no.hyper.real.exposures)
  sp.ovary.parms  <- GetSynSigParamsFromExposures(sp.ovary.exp)


  x.sp.parms <-
    cbind(sp.kidney.parms[ , c("SBS5", "SBS40")],
          sp.ovary.parms[ , "SBS3", drop = FALSE])


  x.sa.parms <-
    cbind(sa.kidney.parms[ , c("BI_COMPOSITE_SBS5_P",
                               "BI_COMPOSITE_SBS40_P")],
          sa.ovary.parms[ , "BI_COMPOSITE_SBS3_P", drop = FALSE])

  sp.abst.info <-
    GenerateSynAbstract(
      parms            = x.sp.parms,
      num.syn.tumors   = num.syn.tumors,
      file.prefix      = NULL, # "sp",
      sample.id.prefix = "SP.Syn.Abst",
      froot            = file.path(top.level.dir, "sp"))

  sa.abst.info <-
    GenerateSynAbstract(
      parms            = x.sa.parms,
      num.syn.tumors   = num.syn.tumors,
      file.prefix      = NULL, #"sa",
      sample.id.prefix = "SA.Syn.Abst",
      froot            = file.path(top.level.dir, "sa"))

  #### Generate and write SignatureAnalyzer "abstract" 3, 5, 40 catalogs

  CreateAndWriteCatalog(
    sa.COMPOSITE.sigs,
    sa.abst.info$syn.exp,
    dir = NULL,
    WriteCatCOMPOSITE,
    overwrite = overwrite,
    my.dir = file.path(top.level.dir, "sa.sa.COMPOSITE"))

  CreateAndWriteCatalog(
    sa.96.sigs,
    sa.abst.info$syn.exp,
    dir = NULL,
    ICAMS::WriteCatalog,
    overwrite = overwrite,
    my.dir = file.path(top.level.dir, "sa.sa.96"))

  # We need to adjust the signature names in the exposures
  # so they match the signature names in \code{sa.COMPOSITE.sigs}.

  tmp.exp <- sp.abst.info$syn.exp
  rownames(tmp.exp) <- rownames(sa.abst.info$syn.exp)

  CreateAndWriteCatalog(
    sa.COMPOSITE.sigs,
    tmp.exp,
    dir = NULL,
    WriteCatCOMPOSITE,
    overwrite = overwrite,
    my.dir = file.path(top.level.dir, "sp.sa.COMPOSITE"))

  CreateAndWriteCatalog(
    sp.sigs,
    sp.abst.info$syn.exp,
    dir = NULL,
    ICAMS::WriteCatalog,
    overwrite = overwrite,
    my.dir = file.path(top.level.dir, "sp.sp"))

  if (!is.null(regress.dir)) {
    diff.result <-
      NewDiff4SynDataSets(top.level.dir, regress.dir, unlink = FALSE)
    if (diff.result[1] != "ok") {
      message("\nThere was a difference, investigate\n",
              paste0(diff.result, "\n"))
    } else {
      message("\nok\n")
    }
  }
}


Create.2.7a.7b.Abstract <- function(seed        = 44,
                                   overwrite   = TRUE,
                                   regress.dir = "data-raw/long.test.regression.data/syn.2.7a.7b.abst/") {

  top.level.dir <- paste0("../syn.3.5.40.abst.", seed)

  set.seed(seed)
  num.syn.tumors <- 1000

  sa.kidney.exp <- GetExpForOneCancerType("Kidney-RCC",
                                          SynSigGen::sa.no.hyper.real.exposures)
  sa.kidney.parms <- GetSynSigParamsFromExposures(sa.kidney.exp)
  sa.ovary.exp <- GetExpForOneCancerType("Ovary-AdenoCA",
                                         SynSigGen::sa.no.hyper.real.exposures)
  sa.ovary.parms  <- GetSynSigParamsFromExposures(sa.ovary.exp)


  sp.kidney.exp <- GetExpForOneCancerType("Kidney-RCC",
                                          SynSigGen::sp.no.hyper.real.exposures)
  sp.kidney.parms <- GetSynSigParamsFromExposures(sp.kidney.exp)
  sp.ovary.exp <- GetExpForOneCancerType("Ovary-AdenoCA",
                                         SynSigGen::sp.no.hyper.real.exposures)
  sp.ovary.parms  <- GetSynSigParamsFromExposures(sp.ovary.exp)


  x.sp.parms <-
    cbind(sp.kidney.parms[ , c("SBS5", "SBS40")],
          sp.ovary.parms[ , "SBS3", drop = FALSE])


  x.sa.parms <-
    cbind(sa.kidney.parms[ , c("BI_COMPOSITE_SBS5_P",
                               "BI_COMPOSITE_SBS40_P")],
          sa.ovary.parms[ , "BI_COMPOSITE_SBS3_P", drop = FALSE])

  sp.abst.info <-
    GenerateSynAbstract(
      parms            = x.sp.parms,
      num.syn.tumors   = num.syn.tumors,
      file.prefix      = NULL, # "sp",
      sample.id.prefix = "SP.Syn.Abst",
      froot            = file.path(top.level.dir, "sp"))

  sa.abst.info <-
    GenerateSynAbstract(
      parms            = x.sa.parms,
      num.syn.tumors   = num.syn.tumors,
      file.prefix      = NULL, #"sa",
      sample.id.prefix = "SA.Syn.Abst",
      froot            = file.path(top.level.dir, "sa"))

  #### Generate and write SignatureAnalyzer "abstract" 3, 5, 40 catalogs

  CreateAndWriteCatalog(
    sa.COMPOSITE.sigs,
    sa.abst.info$syn.exp,
    dir = NULL,
    WriteCatCOMPOSITE,
    overwrite = overwrite,
    my.dir = file.path(top.level.dir, "sa.sa.COMPOSITE"))

  CreateAndWriteCatalog(
    sa.96.sigs,
    sa.abst.info$syn.exp,
    dir = NULL,
    ICAMS::WriteCatalog,
    overwrite = overwrite,
    my.dir = file.path(top.level.dir, "sa.sa.96"))

  # We need to adjust the signature names in the exposures
  # so they match the signature names in \code{sa.COMPOSITE.sigs}.

  tmp.exp <- sp.abst.info$syn.exp
  rownames(tmp.exp) <- rownames(sa.abst.info$syn.exp)

  CreateAndWriteCatalog(
    sa.COMPOSITE.sigs,
    tmp.exp,
    dir = NULL,
    WriteCatCOMPOSITE,
    overwrite = overwrite,
    my.dir = file.path(top.level.dir, "sp.sa.COMPOSITE"))

  CreateAndWriteCatalog(
    sp.sigs,
    sp.abst.info$syn.exp,
    dir = NULL,
    ICAMS::WriteCatalog,
    overwrite = overwrite,
    my.dir = file.path(top.level.dir, "sp.sp"))

  if (!is.null(regress.dir)) {
    diff.result <-
      NewDiff4SynDataSets(top.level.dir, regress.dir, unlink = FALSE)
    if (diff.result[1] != "ok") {
      message("\nThere was a difference, investigate\n",
              paste0(diff.result, "\n"))
    } else {
      message("\nok\n")
    }
  }
}



#' Standardize SignatureAnalyzer signature names.
#'
#' For example, change \code{BI_COMPOSITE_SNV_SBS83_P}
#' to \code{BI_COMPOSITE_SBS83_P}
#'
#' This is necessary because
#' for COMPOSITE signatures we rbind coordinated
#' "SNV", "DNP", and "INDEL" signatures.
#'
#' This is a copy of FixSASigNames in SynSigEval ...
#'
#' @param sig.names Vector of signature names
#'
#' @return Vector of signatures names with "_SNV" removed.
#'
#' @keywords internal

xFixSASigNames <- function(sig.names) {
  return(gsub("_SNV_", "_", sig.names, fixed = TRUE))
}


Defunct.How.We.Mapped.From.SP.to.SA.Sig <-
  function(rcc.sa.syn.exp, ovary.sa.syn.exp) {
  # Find mapping from SBS3, SBS5, and SBS40 to SignatureAnalyzer signatures
  # assigned to these tumor types


  MatchSigs1Direction(
    SynSigGen::sp.sigs[ , "SBS5", drop = F],
    SynSigGen::sa.96.sigs[ , xFixSASigNames(rownames(rcc.sa.syn.exp))])

  MatchSigs1Direction(
    SynSigGen::sp.sigs[ , "SBS40", drop = F],
    SynSigGen::sa.96.sigs[ , xFixSASigNames(rownames(rcc.sa.syn.exp))])

  MatchSigs1Direction(
    SynSigGen::sp.sigs[ , "SBS3", drop = F],
    SynSigGen::sa.96.sigs[ , xFixSASigNames(rownames(ovary.sa.syn.exp))])

  # Both BI..SBS3 and BI..SBS39 are in every ovarian; we select BI..SBS3
}
