#' @title Extract  parameters for one mutational signature profile
#'
#' @param counts  A vector of mutation counts attributed to one signature across
#'                length(counts) samples. TODO(Steve): rename to exposures
#'
#' @param target.size The length of genomic sequence from which the counts
#'                were derived, in megabases. Deprecated, set this to 1.
#'
#' @return
#'   A 3-element vector with names "prevalence", "mean", and "stdev"
#'
#' @importFrom stats sd
#'
#' @keywords internal

SynSigParamsOneSignature <- function(counts, target.size = 1 ) {

  prevalence <-  length(counts[counts >= 1 ]) / length(counts)

  counts.per.mb <- counts[counts >= 1 ] / target.size

  ## generate log10(mut/mb) values for mean and sd
  mean.per.mb <- mean(log10(counts.per.mb))
  sd.per.mb <- sd(log10(counts.per.mb))

  c(prob = prevalence, mean = mean.per.mb, stdev = sd.per.mb)
}


#' @title Empirical estimates of key parameters describing exposures due to signatures.
#'
#' @param exposures A matrix in which each column is a sample and each row is a mutation
#'         signature, with each element being the "exposure",
#'         i.e. mutation count attributed to a
#'         (sample, signature) pair.
#'
#' @param verbose If > 0 cat various messages.
#'
#' @return A data frame with one column for
#' each of a subset of the input signatures
#' and the following rows
#' \enumerate{
#' \item the proportion of tumors with the signature
#' \item mean(log_10(mutations.per.Mb))
#' \item stdev(log_10(mutations.per.Mb))
#' }
#' Signatures not present in
#'  \code{exposures} or present only in a single tumor in
#'  \code{exposures} are removed.
#'
#' @export

GetSynSigParamsFromExposures <- function(exposures, verbose = 0) {
  stopifnot(ncol(exposures) > 0)

  integer.counts <- round(exposures, digits = 0)
  integer.counts <- integer.counts[rowSums(integer.counts) > 0 , ]
  ret1 <- apply(X      = integer.counts,
                MARGIN = 1,
                FUN    = SynSigParamsOneSignature)

  # Some standard deviations can be NA (if there is only one tumor
  # with mutations for that signature). We pretend we did not see
  # these signatures. TODO(Steve): impute from similar signatures.
  if (any(is.na(ret1['stdev', ]))) {
    if (verbose > 0) {
      cat("\nWarning, some signatures present in only one sample, dropping:\n")
      cat(colnames(ret1)[is.na(ret1['stdev', ])], "\n")
    }
  }
  retval <- ret1[,!is.na(ret1['stdev',]) , drop = FALSE]
  if (ncol(retval) == 0) {
    stop("No signatures with usable parameters (> 1 sample with exposure)")
  }
  return(retval)
}

#' @title Write key parameters describing exposures due to a signature to a file.
#'
#' The parameters written are prevalence, mean(log(exposure)), and sd(log(exposure)).
#'
#' @param params The parameters to write.
#'
#' @param file   The path to the file to write.
#'
#' @param append Whether to append to or overwrite \code{file} if it already
#'        exists.
#'
#' @param col.names If NA, add column names.
#'
#' @importFrom utils write.table
#'
#' @export

# Needs to be exported for e.g. Create.pancreas.Rmd

WriteSynSigParams <-
  function(params, file, append = FALSE,
           col.names = ifelse(append, FALSE, NA)) {
    # Suppress warning about writing column names
    # on an append.
    suppressWarnings(
      write.table(x = as.data.frame(params),
                  file = file,
                  sep = ",",
                  col.names = col.names,
                  row.names = TRUE,
                  append = append)
    )
  }

#' @title Create synthetic exposures based given parameters
#'
#' @return A matrix with the rows being each signature and the columns being
#' generated samples. Each entry is the count of mutations due to one
#' signature in one sample.
#'
#' @param sig.params Parameters from \code{\link{GetSynSigParamsFromExposures}}
#'   or another source. Should be
#'   a matrix or data frame with one column for
#'   each signature and the following rows:
#' \describe{
#' \item{prob}{The proportion of tumors with the signature.}
#' \item{mean}{The mean(log_10(number of mutations)).}
#' \item{stdev}{The stdev(log_10(number of mutations)).}
#' }
#'   The rownames need to be the column names of a signature
#'   catalog.
#'
#' @param num.samples Number of samples to generate
#'
#' @param name Prefix for sample identifiers in the simulated dataset
#'
#' @export

GenerateSyntheticExposures <-
  function(sig.params,
           num.samples = 10,
           name = 'synthetic') {

    sigs <- colnames(sig.params)
    stopifnot(!is.null(sigs))
    prev.present <- unlist(sig.params['prob', ]) # Note, get a vector
    sig.burden <- sig.params['mean', , drop = FALSE]
    sig.sd <- sig.params['stdev', , drop = FALSE]

    sig.present <- present.sigs(num.samples, prev.present)

    colnames(sig.present) <- paste(name, seq(1, num.samples), sep = '.')

    # Create a synthetic exposures for each column (sample)
    # in sig.present.
    retval <-
      apply(sig.present,
            2,
            GenerateSynExposureOneSample,
            sigs,
            sig.burden, ## burden is in mutation per megabase
            sig.sd)
    return(retval)
  }


#' Randomly assign present / absent to each of a set of signatures, and
#' keep trying until at least one is present
#'
#' @param prev.present  Vector of prevalences,
#'   each the prevalence of 1 mutational
#'   signature.
#'
#' @param verbose If > 0, cat some possibly informative messages
#'
#' @return a vector of 0s and 1s of length
#' \code{length(prev.present)}, and for which
#' \code{sum(prev.present) > 0}.
#'
#' @importFrom stats rbinom
#'
#' @keywords internal
AssignPresentAbsentOneSample <- function(prev.present, verbose = 0) {
  v <- numeric(length(prev.present))
  while (sum(v) < 1) {
    v <- rbinom(length(prev.present), 1, prev.present)
  }
  if (verbose > 0)
    cat("\nAssignPresentAbsentOneSample returning ", v, "\n\n")
  return(v)
}


#' @title Decide which signatures will be present in
#'  the catalogs of synthetic tumors.
#'
#' @param num.tumors Number of tumors to generate
#'
#' @param prev.present Vector of prevalences,
#'   each the prevalence of 1 mutational
#'   signature. The names are the names of the
#'   signatures.
#'
#' @return A matrix with one row
#' per signature and one column per tumor, with
#' 1 in a cell indicated the presence of a signature
#' and 0 indicating absence.
#'
#' @importFrom stats rbinom
#'
#' @keywords internal

present.sigs <-
  function(num.tumors, prev.present){

    num.process <- length(prev.present)

    present.list <- list()

    for (i in 1:num.process) {
      present.each <- rbinom(num.tumors,
                             1,
                             prob = prev.present[i])
      present.list[[i]] <-  present.each
    }

    present <- do.call(rbind,present.list)

    rownames(present) <- names(prev.present)

    # If the column for one tumor has only 0s,
    # re-sample until there is at least on non-0
    # signature.
    for (tumor in 1:ncol(present)) {
      if (all(present[,tumor] == rep(0, length(nrow(present))))) {
        # present['SBS1',tumor] = 1
        present[ , tumor] <-
          AssignPresentAbsentOneSample(prev.present)
      }
    }
    return(present)
  }


#' @title Using parameters given to generate exposures for synthetic tumors
#'
#' @param tumor Signature presence matrix or exposure matrix for a tumor.
#' It has only one row, and K (# of signatures) columns.
#' Value in each column is the presence flag for a mutational signature:
#' the value can be non-zero(signature is present) or 0(absent).
#' The name of each column should be the name of a signature.
#'
#' @param sig.interest Names of mutational signatures you want to use to
#' generate exposures. It can be all, or part of signatures in colnames(tumor).
#'
#' @param burden.per.sig Mean mutation burden a log10 of the
#' counts of mutations per megabase.
#' It has one row, and K columns. Each column name refers to a mutational signature.
#'
#' @param sd.per.sig standard deviation of mutation burden.
#' It has one row, and K columns. Each column name refers to a mutational signature.
#'
#' @details Determine the intensity of each
#' mutational signature in a tumor, returning the number of mutations
#' using the mean mutation burden per signature and the std dev
#'
#' @importFrom stats rnorm
#'
#' @keywords internal
GenerateSynExposureOneSample <-
  function(tumor,
           sig.interest,
           burden.per.sig,
           sd.per.sig
  ) {

    ## starts with individual tumors, only generate exposures for signatures
    ## with a flag does not equal to 0.
    active.sigs <- base::which(tumor != 0)

    for (sigs in active.sigs) {
      stdev <- sd.per.sig[,sigs]
      burden <- burden.per.sig[,sigs]

      ## if std dev is too big, >= 3, max = 3
      ### consider handling this different. the worry is that the variation
      ##  is too large, the sampled mutation burden will be very high,
      ### which will have a mutation burden that is not biologically possible
      if (stdev >= 3) {
        cat("Very large stdev", stdev, "\n")
        stdev = 3
      }

      ## mutational intensity follows a log normal distibution
      ## use the normal distribution with log-ed values instead
      tumor[sigs] <- 10^(rnorm(1, sd = stdev, mean = burden))
    }

    tumor <- as.matrix(tumor)
    names(tumor) <- sig.interest
    return(tumor)

  }

#' @title Generate synthetic spectra catalogs given
#' signature profiles and synthetic exposures.
#'
#' @param signatures The signature profiles.
#'
#' @param exposures The synthetic exposures.
#'
#' @param sample.id.suffix A string for adding a suffix to
#'  sample ID. For example, if sample.id.suffix is "abc",
#'  then SomeCancerType::s1.33 is changed to
#'  SomeCancerType::s1-abc.33. Actually, this just replaces
#'  the first "." in the sample id with "-" concatenated
#'  to sample.id.suffix. TODO(Steve): probably drop this
#'
#' @return A list of three elements that comprise the
#' synthetic data: \enumerate{
#'  \item \code{ground.truth.catalog}: Spectra catalog for
#'  the software input.
#'  \item \code{ground.truth.signatures}: Signatures active
#'  in \code{ground.truth.catalog}.
#'  \item \code{ground.truth.exposures}: Exposures of \code{ground.truth.signatures}
#'  in \code{ground.truth.catalog}.
#' }
#'
#' @export

CreateSynCatalogs <-
  function(signatures, exposures, sample.id.suffix = NULL) {

  # if (any(colSums(exposures) < 1)) warning("Some exposures < 1")

  exposed.sigs <- rownames(exposures)

  # It is an error if there are signatures in exposures that are not
  # in signatures.
  stopifnot(setequal(setdiff(exposed.sigs, colnames(signatures)), c()))

  # VERY IMPORTANT, the next statement guarantees that
  # the order of signatures in rows of exposures is the same as
  # the order of columns in signatures. In addition,
  # it ensure that signatures contains only signatures
  # that are present in exposures.
  #
  signatures <- signatures[ , exposed.sigs]

  # TODO(Steve): in a future versions, for each
  # synthetic tumor, for each signature, accounting
  # for e mutations in that tumor, sample e mutations
  # from the signature profile.
  catalog <- signatures %*% exposures

  i.cat <- round(catalog, digits = 0)

  if (any(colSums(i.cat) == 0))
    warning("Some tumors with 0 mutations")

  if (!is.null(sample.id.suffix)) {
    newcolnames <-
      gsub(".", paste0("-", sample.id.suffix, "."),
           colnames(i.cat), fixed = TRUE)
    stopifnot(colnames(i.cat == colnames(exposures))) #NEW
    colnames(i.cat) <- newcolnames
    colnames(exposures) <- newcolnames
  }

  ## Convert ground.truth.catalog and ground.truth.signatures
  ## into ICAMS acceptable catalogs before outputting the list
  i.cat <- ICAMS::as.catalog(
    object = i.cat,
    # ref.genome = "hg19",
    # WARNING, there
    abundance = NULL,
      # ICAMS::all.abundance$BSgenome.Hsapiens.1000genomes.hs37d5$genome$`96`,
    region = "genome",
    catalog.type = "counts")

  signatures <- ICAMS::as.catalog(
    object = signatures,
    # ref.genome = "hg19",
    abundance = NULL,
    #  ICAMS::all.abundance$BSgenome.Hsapiens.1000genomes.hs37d5$genome$`96`,
    region = "genome",
    catalog.type = "counts.signature")

  ## Return a list with:
  ## $ground.truth.catalog: Spectra catalog for the software input
  ## $ground.truth.signatures: Signatures active in $ground.truth.catalog
  ## $ground.truth.exposures: Exposures of $ground.truth.signatures in
  ## $ground.truth.catalog.
  return(list(ground.truth.catalog = i.cat,
              ground.truth.signatures = signatures,
              ground.truth.exposures = exposures))
  # TODO(Steve) In future, add noise
  }


#' Merge 2 exposure matrices
#'
#' @param exp1 An exposure matrix
#'
#' @param exp2 An exposure matrix
#'
#' @return The column-wise merge of the two input matrices as
#' with all rownames from either matrix preserved and
#' corresponding entries filled with 0s.
#'
#' @keywords internal
Merge2Exposures <- function(exp1, exp2) {
  # Rows are signatures
  exp.m <- merge(exp1, exp2, by = 0, all = TRUE)
  rownames(exp.m) <- exp.m[ ,1]
  exp.m <- exp.m[ , -1]
  exp.m[is.na(exp.m)] <- 0
  return(as.matrix(exp.m))
}

#' Merge all exposure matrices in a list of matrices
#'
#' @param list.of.exposures A list of exposure matrices
#'
#' @return The column-wise merge of all the input matrices
#' with all rownames from all matrices preserved and
#' corresponding entries filled with 0s.
#'
#' @export

MergeExposures <- function(list.of.exposures) {
  stopifnot(length(list.of.exposures) > 0)
  if (length(list.of.exposures) == 1) {
    return(as.matrix(list.of.exposures[[1]]))
  }
  start <- list.of.exposures[[1]]
  for (i in 2:length(list.of.exposures)) {
    start <- Merge2Exposures(start, list.of.exposures[[i]])
  }
  start2 <- start[order(rownames(start)), ]
  return(as.matrix(start2))
}

#' Create file names in a given directory
#'
#' The directory is provided by the global
#' variable \code{OutDir.dir},
#' which \strong{must} be set by the user. If
#' \code{OutDir.dir} is NULL then just return
#' \code{file.name}.
#'
#' @param file.name The name of the that will be
#' prefixed by \code{OutDir.dir}.
#'
#' @return \code{file.name} prefixed by \code{OutDir.dir}.
#'
#' @export
# Paste of OutDir.dir / file.name (or just file.name)
OutDir <- function(file.name) {
  warning("Use of function OutDir is deprecated")
  if (!exists("OutDir.dir")) return(file.name)
  if (is.null(OutDir.dir)) return(file.name)
  tmp <- OutDir.dir
  n <- nchar(tmp)
  if (substr(tmp, n, n) != "/")
    tmp <- paste0(tmp, "/")
  if (!dir.exists(tmp)) {
    dir.create(tmp)
  }
  return(paste0(tmp, file.name))
}


#' Generate synthetic exposures from abstract parameters.
#'
#' Checkpoints the parameters and the synthetic
#' exposures to files. It also checks that the parameters
#' inferred from the synthetic data approximate those
#' inferred from \code{real.exp}.
#'
#' @param parms The actual exposures upon which to base
#' the parameters and synthetic exposures.
#'
#' @param num.syn.tumors Generate this number of synthetic tumors.
#'
#' @param file.prefix Prepend this to output filenames
#'  to indicate the organization of the data.
#'
#' @param sample.id.prefix Prefix for sample identifiers for the
#' synthetic samples.
#'
#' @return A list with elements:
#' \enumerate{
#'  \item \code{parms} The parameters inferred from \code{real.exp}.
#'  \item \code{syn.exp} The synthetic exposures generated from \code{parms}.
#' }
#'
#' @keywords internal

GenerateSynAbstract <-
  function(parms, num.syn.tumors, file.prefix = NULL, sample.id.prefix, froot = NULL) {
    stopifnot(!is.null(parms))

    if (is.null(froot)) {
      froot <- OutDir(file.prefix)
    }

    parm.file <- paste0(froot, ".parms.csv")
    cat("# Original paramaters\n", file = parm.file)
    suppressWarnings( # Suppress warning on column names on append
      WriteSynSigParams(parms, parm.file, append = TRUE,
                      col.names = NA))

    syn.exp <-
      GenerateSyntheticExposures(parms, num.syn.tumors, sample.id.prefix)

    WriteExposure(syn.exp, paste0(froot, ".generic.syn.exposure.csv"))

    # Sanity check; we regenerate the parameters from the synthetic exposures.
    check.params <- GetSynSigParamsFromExposures(syn.exp)

    # check.params should be similar to parms
    cat("# Parameters derived from synthetic exposures\n",
        file = parm.file, append = TRUE)
    suppressWarnings(
      WriteSynSigParams(check.params, parm.file, append = TRUE))

    missing.sig.names <- setdiff(colnames(parms), colnames(check.params))
    if (length(missing.sig.names) > 0) {
      cat("# Some signatures not represented in the synthetic data:\n",
          file = parm.file, append =  TRUE)
      cat("#", missing.sig.names, "\n", file = parm.file, append = TRUE)
      check.param2 <- matrix(NA, nrow = dim(parms)[1], ncol = dim(parms)[2])
      dimnames(check.param2) <- dimnames(parms)
      check.param2[ , colnames(check.params)] <- check.params
      check.params <- check.param2
    }

    cat("# Difference between original parameters and parameters",
        "derived from synthetic exposures\n",
        file = parm.file, append = TRUE)
    WriteSynSigParams(parms - check.params, parm.file,
                      append = TRUE)
    cat("# The difference should be small\n",
        file = parm.file, append = TRUE)

    return(list(parms = parms, syn.exp = syn.exp, check.parms = check.params))
  }

#' Generate synthetic exposures from real exposures.
#'
#' Checkpoints the parameters and the synthetic
#' exposures to files. It also checks that the parameters
#' inferred from the synthetic data approximate those
#' inferred from \code{real.exp}.
#'
#' @param real.exp The actual (real) exposures upon which to base
#' the parameters and synthetic exposures.
#'
#' @param num.syn.tumors Generate this number of synthetic tumors.
#'
#' @param file.prefix Prepend this to output filenames
#'  to indicate the organization of the data.
#'
#' @param sample.id.prefix Prefix for sample identifiers for the
#' synthetic samples.
#'
#' @param top.level.dir Directory in which to create several files.
#'   This directory must already exist.
#'
#' @return A list with elements:
#' \enumerate{
#'  \item \code{parms} The parameters inferred from \code{real.exp}.
#'  \item \code{syn.exp} The synthetic exposures generated from \code{parms}.
#' }
#'
#' @export
GenerateSynFromReal <- function(real.exp,
                                num.syn.tumors,
                                file.prefix,
                                sample.id.prefix,
                                top.level.dir = NULL) {

  stopifnot(ncol(real.exp) > 0)
  parms <- GetSynSigParamsFromExposures(real.exp)

  if (is.null(top.level.dir)) {
    new.file.prefix <- OutDir(file.prefix)
    warning("Calling GenerateSynFromReal witout top.level.dir is deprecated")
  } else {
    new.file.prefix <- file.path(top.level.dir, file.prefix)
  }

  WriteExposure(real.exp,
                paste0(new.file.prefix, ".real.input.exposure.csv"))

  return(
    GenerateSynAbstract(
      parms = parms,
      num.syn.tumors = num.syn.tumors,
      file.prefix = file.prefix,
      sample.id.prefix = sample.id.prefix,
      froot = file.path(top.level.dir, file.prefix)))

  #   function(parms, num.syn.tumors, file.prefix, sample.id.prefix, froot = NULL)
}

#' Create and write a mutational spectra catalog
#'
#' @export
#'
#' @param sigs Signatures to use.
#'
#' @param exp (Synthetic) exposures.
#'
#' @param dir Deprecated, maintained only to avoid
#' breaking old code. A subdirectory based on
#' the deprecated global variable \code{\link{OutDir}}.
#'
#' @param write.cat.fn Function to write catalogs \strong{or}
#' spectra to files.
#'
#' @param extra.file.suffix Extra string to put before ".csv".
#'
#' @param overwrite If TRUE, overwrite existing directory; useful for
#' debugging / testing.
#'
#' @param my.dir The directory in which to write the catalog
#'  and several additional files.
#'
#' @return Invisibly, the generated catalog.
#'
#' @details Create a file with the catalog \code{syn.data.csv}
#'  and writes \code{sigs} to \code{input.sigs.csv}.
#'
CreateAndWriteCatalog <-
  function(sigs,
           exp,
           dir               = NULL,
           write.cat.fn      = ICAMS::WriteCatalog,
           extra.file.suffix = "",
           overwrite         = FALSE,
           my.dir            = NULL) {
    info <- CreateSynCatalogs(sigs, exp)

    if (is.null(my.dir)) {
      warning("In CreateAndWriteCatalog my.dir is NULL, using deprecated work-around\n",
              "Do not use the argument dir, use out.dir")
      if (exists("OutDir") && exists("OutDir.dir")) {
        my.dir <-  OutDir(dir)
      } else {
        warning("In CreateAndWriteCatalog my.dir is NULL and OutDir() does not exist\n",
                "setting my.dir <- dir")
        my.dir <- dir
      }
    }

    if (dir.exists(my.dir)) {
      if (!overwrite) stop("\nDirectory ", my.dir, " exists\n")
    } else {
      dir.create(my.dir)
    }

    if (extra.file.suffix == "") {
      suffix <- ".csv"
    } else {
      suffix = paste0(".", extra.file.suffix, ".csv")
    }
    write.cat.fn(info$ground.truth.signatures,
                  paste0(my.dir, "/ground.truth.syn.sigs", suffix))

    zero.mutation <- base::which(colSums(info$ground.truth.catalog) == 0)

    if (length(zero.mutation) > 0) {
      warning("Tumors with no mutation:\n\n",
              colnames(info$ground.truth.catalog)[zero.mutation],
              "in", my.dir)
    }
    write.cat.fn(info$ground.truth.catalog,
                 paste0(my.dir, "/ground.truth.syn.catalog", suffix))
    WriteExposure(info$ground.truth.exposures,
                  paste0(my.dir, "/ground.truth.syn.exposures", suffix))
    invisible(info$ground.truth.catalog)
  }

#' @keywords internal
MustCreateDir <- function(dir, overwrite = FALSE) {
  if (dir.exists(dir)) {
    if (!overwrite) stop(dir, " exists and overwrite is FALSE")
    return(NULL)
  }
  if (!dir.create(dir, recursive = TRUE)) {
    stop("Unable to create dir ", dir )
  }
  return(NULL)
}

#' Create and write a mutational spectra catalog
#'
#' @export
#'
#' @param sigs Signatures to use.
#'
#' @param exp (Synthetic) exposures.
#'
#' @param dir Directory in which to put the signatures;
#' NOTE: this will be a subdirectory based on \code{\link{OutDir}}.
#'
#' @param extra.file.suffix Extra string to put before ".csv".
#'
#' @param overwrite If TRUE, overwrite existing directory; useful for
#' debugging / testing.
#'
#' @return Invisibly, the generated catalog.
#'
#' @details Create a file with the catalog \code{syn.data.csv}
#'  and writes \code{sigs} to \code{input.sigs.csv}.
#'
NewCreateAndWriteCatalog <-
  function(sigs, exp, dir, extra.file.suffix = "",
           overwrite = FALSE) {
    info <- CreateSynCatalogs(sigs, exp)

    if (dir.exists(dir)) {
      if (!overwrite) stop("\nDirectory ", dir, " exists\n")
    } else {
      MustCreateDir(dir)
    }

    if (extra.file.suffix == "") {
      suffix <- ".csv"
    } else {
      suffix = paste0(".", extra.file.suffix, ".csv")
    }
    ICAMS::WriteCatalog(info$ground.truth.signatures,
                        paste0(dir, "/ground.truth.syn.sigs", suffix))

    zero.mutation <- base::which(colSums(info$ground.truth.catalog) == 0)

    if (length(zero.mutation) > 0) {
      warning("Tumors with no mutation:\n\n",
              colnames(info$ground.truth.catalog)[zero.mutation],
              "in", dir)
    }
    ICAMS::WriteCatalog(info$ground.truth.catalog,
                        paste0(dir, "/ground.truth.syn.catalog", suffix))
    WriteExposure(info$ground.truth.exposures,
                  paste0(dir, "/ground.truth.syn.exposures", suffix))
    invisible(info$ground.truth.catalog)
  }


#' Exposures and spectra with Poisson or negative binomial noise in exposures.
#'
#' @param input.exposure The exposures to which to add noise; a numeric matrix
#'    or data frame in which the rows are signatures and the columns are
#'    samples. Each cell indicates the number of mutations due to a particular
#'    signature in a particular sample.
#'
#' @param signatures The signatures in the exposure; the column names
#'     of \code{signatures} have to include all row names in
#'     \code{input.exposure}; can be an \code{\link[ICAMS]{ICAMS}}
#'     catalog or a numerical matrix or data frame.
#'
#' @param n.binom.size If non \code{NULL}, use negative binomial noise
#'     with this size parameter; see \code{\link[stats]{NegBinomial}}.
#'
#' @return A list with the elements \describe{
#' \item{expsoures}{The numbers of mutations due to each signature
#'    after adding noise}
#' \item{spectra}{The spectra based on the noisy signature exposures.}
#' }
#'
#' @importFrom stats rpois rnbinom
#'
#' @export

AddNoise <- function(input.exposure, signatures, n.binom.size = NULL) {

  if (is.data.frame(input.exposure)) {
    input.exposure <- as.matrix(input.exposure)
  }
  stopifnot(is.numeric(input.exposure))
  exposures <- input.exposure
  exposures[ , ] <- NA

  spectra <- matrix(0, ncol = ncol(input.exposure), nrow = nrow(signatures))
  rownames(spectra) <- rownames(signatures)
  colnames(spectra) <- colnames(input.exposure)

  for (sig in rownames(input.exposure)) {

    partial.spec <- signatures[ , sig] %*% input.exposure[sig, , drop = FALSE]
    # Resample (add noise) to the "partial spectra" due to sig
    if (is.null(n.binom.size)) {
      noised.vec <-
        rpois(n = length(partial.spec), lambda = partial.spec)
    } else {
      noised.vec <-
        rnbinom(n = length(partial.spec), size = n.binom.size, mu = partial.spec)
    }
    # Turn the vector back into a matrix
    noisy.partial.spec <- matrix(noised.vec, nrow = nrow(partial.spec))
    exposures[sig, ] <- colSums(noisy.partial.spec)
    spectra <- spectra + noisy.partial.spec

  }

  return(list(exposures = exposures, spectra = spectra))

}
