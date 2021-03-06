### A script to generate synthetic exposures for TWO mutational signatures in
### the log space using parameters (prob = 1, mean.log, main.stdev.log) for TWO
### mutational signatures of interest Note: This generator generates samples
### whose exposures of two signatures are PARTIALLY CORRELATED, This is
### realized by letting the mean correlated.signature exposure to be the
### exact value of main.signature exposure for each data point.

#' @title Generate correlated exposures for one tumor
#'
#' Function to generate exposure of two correlated signatures
#' (Example: SBS1 and SBS5) for ONE synthetic tumor.
#'
#' @param tumor.name Name of synthetic tumor you want to generate.
#' (default: \code{TwoCorreSigsGen::1})
#'
#' @param main.signature Name of a signature with smaller variance
#' in the log10 space. (default: "SBS5")
#'
#' @param correlated.signature Name of a signature with larger variance
#' in the log10 space. (default: "SBS1")
#'
#' @param main.mean.log Mean of log10(mutation burden of \code{main.signature})
#'
#' @param main.stdev.log Standard deviation of log10(mutation burden
#' of \code{main.signature})
#'
#' @param correlated.stdev.log Contribute to part of the standard deviation of
#' log10(mutation burden of \code{correlated.signature}). In this script, the s.d. of
#' log10(mutation burden of \code{correlated.signature})
#' = \code{main.stdev.log} + \code{correlated.stdev.log}
#'
#' @param slope.linear Average ratio of mutation burden of \code{correlated.signature}
#' over mutation burden of \code{main.signature}
#'
#' @param main.signature.lower.thres Minimum mutation burden (number of mutations)
#' induced by \code{main.signature} in each tumor.
#'
#' @param correlated.signature.lower.thres Minimum mutation burden (number of mutations)
#' induced by \code{correlated.signature} in each tumor.
#'
#' @param min.main.to.correlated.ratio.linear Minimum ratio of
#' \code{main.signature} over mutation burden of
#' \code{correlated.signature} in each tumor.
#'
#' @param max.main.to.correlated.ratio.linear Maximum ratio of
#' \code{main.signature} over mutation burden of
#' \code{correlated.signature} in each tumor.
#'
#' @section Warning Warning 1: forcing dataset to have high Pearson's R^2,
#' while inputting large main.stdev.log value will cause
#' the program to reject unqualified datasets and regenerate
#' new datasets for HOURS OR EVEN DAYS.
#'
#' Warning 2: If users intend to replicate results of SBS1-SBS5
#' paper, please specify parameters as indicated in the
#' Supplementary Information of the paper.
#'
#' @keywords internal
#'
GenSBS1SBS5ExposureOneTumor <- function(
  tumor.name = "TwoCorreSigsGen::1",
  main.signature = "SBS5",
  correlated.signature = "SBS1",
  main.mean.log = 2.5,
  main.stdev.log = 0.25,
  correlated.stdev.log = 0.25,
  slope.linear = 1,
  main.signature.lower.thres = 50,
  correlated.signature.lower.thres = 30,
  min.main.to.correlated.ratio.linear = 1/3,  ## The lower ratio for count(SBS5) / count(SBS1) in LINEAR SPACE!
  max.main.to.correlated.ratio.linear = Inf)  ## The higher ratio for count(SBS5) / count(SBS1) in LINEAR SPACE!
{
  #### Matrices storing exposures for two signatures
  exposure.counts <- matrix(data = NA,              ## Referring to test.exposure.counts,
                            nrow = 2,               ## Each row refers to a signature, (Thus nrow = 2)
                            ncol = 1)               ## And each column refers to a sample. (Thus ncol = 1)


  #### Add names for signatures of interest and tumors to be generated
  rownames(exposure.counts) <- c(main.signature,correlated.signature)
  colnames(exposure.counts) <- tumor.name

  ## Repeat generating exposure for both,
  ## Until the count is greater then lower.thres AND
  ## The count(main)/count(correlated) ratio is between two thresholds
  repeat{

    ## Repeat generating the exposure for main.signature,
    ## until the count is greater than main.signature.lower.thres
    repeat{

      ## First construct a signature presence matrix for main.signature
      sig.presence <- matrix(data = c(1), nrow = 1, ncol = 1)
      rownames(sig.presence) = tumor.name
      colnames(sig.presence) = main.signature

      ## Construct mean.burden.per.sig matrix for main.signature
      mean.per.sig <- matrix(main.mean.log,1,1)
      rownames(mean.per.sig) = tumor.name
      colnames(mean.per.sig) = main.signature

      ## Construct sd.per.sig for main.signature
      sd.per.sig <- matrix(main.stdev.log,1,1)
      rownames(sd.per.sig) = tumor.name
      colnames(sd.per.sig) = main.signature

      ## Using GenerateSynExposureOneSample() to generate exposure
      ## of main.signature for a tumor.
      exposure.counts[main.signature,1] <-
        GenerateSynExposureOneSample(tumor = sig.presence,
                                     sig.interest = main.signature,
                                     burden.per.sig = mean.per.sig,
                                     sd.per.sig = sd.per.sig)
      main.signature.count.log <- log10(exposure.counts[main.signature,1])

      if(exposure.counts[main.signature,1] >= main.signature.lower.thres )
        break  ## Keep regenerating the data, until current generated exposure is greater than threshold
    }

    ## Repeat generating the exposure for correlated.signature,
    ## until the count is greater than correlated.signature.lower.thres
    repeat{

      ## Next construct a signature presence matrix for correlated.signature
      sig.presence <- matrix(data = c(1), nrow = 1, ncol = 1)
      rownames(sig.presence) = tumor.name
      colnames(sig.presence) = correlated.signature

      ## Construct mean.burden.per.sig matrix for correlated.signature
      ## To enable strong correlation between exposures of two signatures,
      ## For each tumor, the log-mean burden of correlated signature
      ## equals to log value of main signature in this tumor + log10(slope)
      mean.per.sig <- matrix(main.signature.count.log+log10(slope.linear),1,1)
      rownames(mean.per.sig) = tumor.name
      colnames(mean.per.sig) = correlated.signature

      ## Construct sd.per.sig for correlated.signature
      sd.per.sig <- matrix(correlated.stdev.log,1,1)
      rownames(sd.per.sig) = tumor.name
      colnames(sd.per.sig) = correlated.signature

      ## Using GenerateSynExposureOneSample() to generate exposure
      ## of correlated.signature for a tumor.
      exposure.counts[correlated.signature,1] <-
        GenerateSynExposureOneSample(tumor = sig.presence,
                                     sig.interest = correlated.signature,
                                     burden.per.sig = mean.per.sig,
                                     sd.per.sig = sd.per.sig)


      if( exposure.counts[correlated.signature,1] >= correlated.signature.lower.thres )
        break  ## Keep regenerating the data, until all generated exposures are greater than threshold
    }

    main.to.correlated.ratio.linear <- exposure.counts[main.signature,1] / exposure.counts[correlated.signature,1]

    if( main.to.correlated.ratio.linear >= min.main.to.correlated.ratio.linear & main.to.correlated.ratio.linear <= max.main.to.correlated.ratio.linear )
      break

  }

  return(exposure.counts)

}


## Generic function for generating signature exposure matrix.
## sample.number specify the number of tumors in this matrix
## Parameters control tumors' exposure amount to the main.signature(SBS5) and correlated signature(SBS1)
## This function is defined for tuning the parameters, and reshape the datasets more flexibily
## Signature exposure matrix can be served as "ground-truth" or "benchmark" attribution,
## and thus can be used to compare with attribution results from different software tools.

#' @title Generate correlated exposures for multiple tumors
#'
#' Wrapper function around \code{\link{GenSBS1SBS5ExposureOneTumor}}:
#' A function to generate exposure of two correlated signatures
#' (Example: SBS1 and SBS5) for \code{sample.number} (e.g. 500) synthetic tumors.
#'
#' NOTE: \code{pearson.r.2.lower.thres} and \code{pearson.r2.higher.thres}
#' are used to constraint the Pearson's R^2 of mutation burdens of two signatures
#' in multiple tumors.
#'
#'
#' @param main.signature Name of a signature with smaller variance
#' in the log10 space. (Default: "SBS5")
#'
#' @param correlated.signature Name of a signature with larger variance
#' in the log10 space. (Default: "SBS1")
#'
#' @param sample.number Number of tumors whose mutation burdens
#' will be generated. (Default: 500)
#'
#' @param name.prefix Prefix of tumor name.
#' (Default: \code{"TwoCorreSigsGen"})
#' By default, the name of tumors to be created will be:
#' TwoCorreSigGen::1, TwoCorreSigGen::2, TwoCorreSigGen::3...
#'
#' @param main.mean.log Mean of log10(mutation burden of \code{main.signature})
#'
#' @param main.stdev.log Standard deviation of log10(mutation burden
#' of \code{main.signature})
#'
#' @param correlated.stdev.log Contribute to part of the standard deviation of
#' log10(mutation burden of correlated.signature). In this script, the s.d. of
#' log10(mutation burden of correlated.signature)
#' = main.stdev.log + correlated.stdev.log
#'
#' @param slope.linear Average ratio of mutation burden of \code{correlated.signature}
#' over mutation burden of \code{main.signature}
#'
#' @param main.signature.lower.thres Minimum mutation burden
#' (number of mutations) induced by \code{main.signature} in each tumor.
#'
#' @param correlated.signature.lower.thres Minimum mutation burden
#' (number of mutations) induced by \code{correlated.signature} in each tumor.
#'
#' @param pearson.r.2.lower.thres Minimum Pearson's R^2 of mutation burdens
#' of two signatures in \code{sample.number} tumors.
#'
#' @param pearson.r.2.higher.thres Maximum Pearson's R^2 of mutation burdens
#' of two signatures in \code{sample.number} tumors.
#'
#' @param min.main.to.correlated.ratio.linear Minimum ratio of
#' \code{main.signature} over mutation burden of
#' \code{correlated.signature} in each tumor.
#'
#' @param max.main.to.correlated.ratio.linear Maximum ratio of
#' \code{main.signature} over mutation burden of
#' \code{correlated.signature} in each tumor.
#'
#' @importFrom stats cor
#'
#' @export
#'
GenSBS1SBS5Exposure <- function(
  main.signature = "SBS5",
  correlated.signature = "SBS1",
  sample.number = 500,
  name.prefix = "TwoCorreSigsGen",
  main.mean.log = 2.5,
  main.stdev.log = 0.25,
  correlated.stdev.log = 0.25,
  slope.linear = 1,
  main.signature.lower.thres = 50,
  correlated.signature.lower.thres = 30,
  pearson.r.2.lower.thres = 0.1,
  pearson.r.2.higher.thres = 1.0,
  min.main.to.correlated.ratio.linear = 1/3,    ## The lower ratio for count(SBS5) / count(SBS1) in LINEAR SPACE!, minimum value is 0
  max.main.to.correlated.ratio.linear = Inf) ## The higher ratio for count(SBS5) / count(SBS1) in LINEAR SPACE! maximum value is Inf
{
  #### Matrices storing exposures for two signatures
  #### For a specific tumor sample, exposure.mb records exposure number of a signature per unit megabase
  #### exposure.counts records exposure number of a signature throughout whole genome/exome
  #### We use exposure.counts to test the performance of differnet tools.
  exposure.counts <- matrix(data = NA,              ## Referring to test.exposure.counts,
                            nrow = 2,               ## Each row refers to a signature, (Thus nrow = 2)
                            ncol = sample.number)   ## And each column refers to a sample. (Thus ncol = sample.number)


  #### Add names for signatures of interest and tumors to be generated
  rownames(exposure.counts) <- c(main.signature,correlated.signature)
  colnames(exposure.counts) <- paste(name.prefix,seq(1,sample.number),sep = "::") ## Designate names for synthetic tumors

  repeat{
    #### Repeat generating exposures for both signatures
    #### Note: For each individual tumor, we need to discard the generated exposure,
    #### if the count(main.signature) < main.signature.lower.thres
    for( ii in 1:sample.number ){ ## Generate a single tumor at a time in order to let the discarding process faster
      temp.exposure.counts <-
        GenSBS1SBS5ExposureOneTumor(tumor.name = paste0(name.prefix,"::",ii),
                                    main.signature = main.signature,
                                    correlated.signature = correlated.signature,
                                    main.mean.log = main.mean.log,
                                    main.stdev.log = main.stdev.log,
                                    correlated.stdev.log = correlated.stdev.log,
                                    slope.linear = slope.linear,
                                    main.signature.lower.thres = main.signature.lower.thres,
                                    correlated.signature.lower.thres = correlated.signature.lower.thres,
                                    min.main.to.correlated.ratio.linear = min.main.to.correlated.ratio.linear,   ## The lower ratio for count(SBS5) / count(SBS1) in LINEAR SPACE!
                                    max.main.to.correlated.ratio.linear = max.main.to.correlated.ratio.linear)

      exposure.counts[,ii] <- temp.exposure.counts
    }

    a <- exposure.counts[main.signature,]
    b <- exposure.counts[correlated.signature,]
    ## If Pearson's R^2 > pearson.r.2.lower.thres in log-transformed dataset, then accept this dataset
    if( pearson.r.2.higher.thres >= cor(log10(a),log10(b)) ^ 2  & cor(log10(a),log10(b)) ^ 2 > pearson.r.2.lower.thres )
      break
  }

  return(exposure.counts)
}


#' @title Plot scatter plot for correlation between two vectors.
#'
#'
#' \code{PlotCorrelationScatterplot} is a wrapper around \code{graphics::plot()},
#' and a function to plot the correlation between two vectors,
#' \code{x} and \code{y}. These vectors are expected to be
#' exposures of two signatures.
#'
#' It will draw a scatterplot, and it will also print information
#' onto the plot, including correlation between \code{x} and \code{y},
#' mean and stdev of \code{x} and \code{y}, etc.
#'
#'
#' @param x vector of exposures of \code{main.signature}
#' (SBS5 in the paper).
#' The exposures of \code{main.siganture} will be aligned onto x axis.
#'
#' @param y vector of exposures of \code{correlated.signature}
#' (SBS1 in the paper).
#' The exposures of \code{correlated.siganture} will be aligned onto y axis.
#'
#' @param xlab Label below x axis.
#'
#' @param ylab Label below y axis.
#'
#' @param main Title on the scatterplot.
#' Default: NULL
#'
#' @param optional.remarks Remarks added below the title.
#'
#' @param ... Other parameters provided to the function \code{graphics::plot()}.
#'
#' @importFrom stats cor var
#' @importFrom graphics mtext plot
#'
#'
#' @export
#'

PlotCorrelationScatterplot <- function(
  x,y,
  xlab = NULL,
  ylab = NULL,
  main = NULL,
  optional.remarks = "",
  ...) {
  ## Check whether x and y are vectors of same length
  stopifnot(is.vector(x) & is.vector(y))
  stopifnot(length(x) == length(y))

  plot(x = x, y = y,
       xlab = xlab,
       ylab = ylab,
       main = main,
       ... = ...)

  mtext(paste("Pearson R^2 = ", round(cor(x,y)^2,3), "; ",
              "x.stdev = ", round(sqrt(var(x)),3), "; ",
              "y.stdev = ", round(sqrt(var(y)),3), "; ",
              sep = ""),
        cex = 0.8) ## By default, the mtext() function adds text at the top margin
  mtext(paste("x.mean = ", round(mean(x),3), "; ",
              "y.mean = ", round(mean(y),3), "; ",
              "Number of data points: ", length(x),
              sep = ""),
        line = 1,
        cex = 0.8) ## By default, the mtext() function adds text at the top margin
  mtext(optional.remarks, ## Adds a line of optional.remarks above the previous text line.
        line = 2,
        cex = 0.8) ## By default, the mtext() function adds text at the top margin
}


#' @title Plot scatter plot for correlation between exposures of two signatures
#'
#' Plot scatter plot for correlation between exposures of two signatures,
#' SBS1 and SBS5 in this study.
#'
#' \code{PlotCorrelationScatterplotForExposures}
#' is a wrapper around \code{\link{PlotCorrelationScatterplot}}.
#' It lets \code{exposure.counts <-} the exposure matrix,
#' and will draw a scatterplot for exposures of two signatures.
#'
#'
#' @param main Title on the scatterplot.
#' Default: NULL
#'
#' @param pdf.filename Name of the PDF to contain the scatterplots.
#'
#' @param main.signature Name of a signature with smaller variance
#' in the log10 space. (Default: "SBS5")
#'
#' @param correlated.signature Name of a signature with larger variance
#' in the log10 space. (Default: "SBS1")
#'
#' @param slope.linear Average ratio of mutation burden of \code{correlated.signature}
#' over mutation burden of \code{main.signature}
#'
#' @param exposure.counts Data.frame or matrix storing exposures of two signatures.
#' The exposure.counts object is usually obtained from \code{SynSig::ReadExposure()}.
#'
#' @param xlim,ylim numeric vectors of length 2,
#' giving the x and y coordinates ranges.
#' Default: c(0,4)
#'
#' @param ... Other parameters provided to the function \code{graphics::plot()}.
#'
#' @importFrom grDevices dev.off pdf
#'
#' @export
#'
PlotCorrelationScatterplotForExposures <-
  function(pdf.filename,
           main.signature = "SBS5",
           correlated.signature = "SBS1",
           slope.linear,
           exposure.counts,
           xlim=c(0,4),
           ylim=c(0,4),
           ...) {
    pdf(pdf.filename)

    PlotCorrelationScatterplot(
      x = log(exposure.counts[main.signature,], base = 10),
      y = log(exposure.counts[correlated.signature,], base = 10),
      xlab = paste("log10( ",main.signature," )",sep = ""),
      ylab = paste("log10( ",correlated.signature," )",sep = ""),
      main = "",
      optional.remarks = paste("slope.linear = ", slope.linear,"; ",
                               sep = ""),
      xlim = xlim,
      ylim = ylim,
      ... = ...)

    dev.off()
  }

##########################################################################################
######## Generate Dataset in which the exposure of two signatures are correlated.
######## Generate the exposure for SBS5 (main.signature) first,
######## Then generate SBS1 (correlated.signature) exposure correlated to the SBS5 exposure.
##########################################################################################

##########################################################################################

#' Wrapper function for generating SBS1-SBS5-correlated Synthetic data
#'
#' This function will use SigProfiler-SBS96 mutational signatures
#' to generate imaginary tumor spectra with mutation burdens only
#' from SBS1 and SBS5, and mutation burdens of both signatures
#' are highly correlated.
#'
#' If you want to customize Pearson R^2 of the dataset,
#' you need to change the standard deviations of two signatures.
#' i.e., main.stdev.log and correlated.stdev.log.
#'
#' This function will generate files listed below:
#'
#' ground.truth.syn.catalog.csv: Generated tumor spectra in
#' ICAMS SBS96 CSV format.
#'
#' ground.truth.syn.exposures.csv: Mutation burdens of SBS1 and
#' SBS5 in generated tumor spectra in ICAMS CSV format.
#'
#' ground.truth.syn.sigs.csv: Ground-truth SBS1 and SBS5
#' signatures in ICAMS SBS96 CSV format.
#'
#' parameters.txt: Parameters used to generate the exposures
#' and tumor spectra.
#'
#' scatterplot.pdf: scatterplot illustrating correlation of
#' exposures of two signatures in generated spectra
#'
#' seedInUse.txt, RNGInUse.txt: seed and Random Number Generator
#' used in generation. (For better reproducibility)
#'
#' sessionInfo.txt: information related to R versions, platforms,
#' loaded or imported packages, etc. (For better reproducibility)
#'
#'
#' @param dir.name Folder to place the generated tumor spectra
#' and other output files.
#' Default: ./S.0.5.Rsq.0.3
#'
#' @param dataset.name The dataset.name encodes the parameters for
#' the synthetic data, but this is just a convention.
#' If NULL, it will be changed to the last part of the \code{dir.name}
#' (Default: NULL)
#'
#' @param overwrite Whether to overwrite
#' (Default: FALSE)
#'
#' @param seed The seed number used to initialize pseudo-random number
#' generator (RNG). This makes the generation of the correlated
#' datasets repeatable. (Default: 1)
#'
#' @param parameter.df a named 1*14 data.frame containing the following items:
#'
#' \enumerate{
#'
#' \item \code{main.signature} The name of the main signature whose exposure
#' can vary freely. (Default: SBS5)
#'
#' \item \code{correlated.signature} The name of the correlated signature
#' whose exposure is influenced by and co-varies with the exposure
#' of main.signature. In this study, it defaults as "SBS1".
#'
#' \item \code{name.prefix} Default: \code{TwoCorreSigsGen}
#'
#' \item \code{sample.number} The number of synthetic tumors you want to generate.
#' Default: 500
#'
#' \item \code{main.mean.log} The mean of log(count(SBS5),base = 10)
#' Default: 2.5
#'
#' \item \code{main.stdev.log} The standard deviation of log(count(SBS5),base = 10)
#' Default: 0.3
#'
#' \item \code{correlated.stdev.log} The ADDED standard deviation of log(count(SBS1),base = 10).
#' This parameter is ADDED stdev because based on the mechanism to generate the count,
#' log10(count(SBS1)) inherently has a stdev = slope * main.stdev.log
#' Default: 0.4
#'
#' \item \code{slope.linear} The ratio for: (Correlated exposure) / (Main exposure) IN LINEAR SPACE!
#' Default: 0.5
#'
#' \item \code{main.signature.lower.thres} This program will force the exposure count of
#' main.signature to be greater than this threshold.
#' Default: 100
#'
#' \item \code{correlated.signature.lower.thres} This program will force the exposure count of
#' correlated.signature to be greater than this threshold.
#' Default: 1
#'
#' \item \code{pearson.r.2.lower.thres} Lower boundary of Pearson's R^2
#' (Default: 0.29)
#'
#' \item \code{pearson.r.2.higher.thres} Upper boundary of Pearson's R^2
#' (Default: 0.31)
#'
#' \item \code{min.main.to.correlated.ratio.linear} The lower ratio for count(SBS5) / count(SBS1)
#' in LINEAR SPACE! (Default: 1/3)
#'
#' \item \code{max.main.to.correlated.ratio.linear} The upper ratio for count(SBS5) / count(SBS1)
#' in LINEAR SPACE! (Default: Inf)
#'
#' }
#'
#' @param add.info Whether to generate additional information.
#'
#' @param verbose If \code{TRUE} cat progress messages.
#' You should set it to \code{FALSE} when you want to make a \code{diff}
#' using \code{\link{CreateSBS1SBS5CorrelatedSyntheticData}}
#' (i.e. parameter \code{regressdir} is not \code{NULL}). This is because
#' Additional information may differ on different OS or R sessions,
#' thus may prevent the dataset from passing the \code{NewDiff4SynDatasets} check.
#' (Default: TRUE)
#'
#' \strong{Warning} \cr
#' Exposure generation function will repeat generating exposure counts
#' using mean and stdev parameters, until the dataset has a Pearson's R^2
#' which falls between two boundaries of Pearson's R^2.
#' Below are a group of parameters which have been tested successfully.
#' If you intend to lower the Pearson's R^2, do remember to increase
#' the main.stdev.log and correlated.stdev.log.
#' Otherwise, the exposure generation will keep generating and discarding datasets!
#'
#' @importFrom ICAMS WriteCatalog
#' @importFrom utils capture.output sessionInfo
#'
#' @export
#'
CreateSBS1SBS5CorrelatedSyntheticDataOneDataset <-
  function(dir.name = "./S.0.5.Rsq.0.3",
           dataset.name = NULL,
           overwrite = FALSE,
           seed = 1,
           parameter.df = SynSigGen::SBS1SBS5parameter["S.0.5.Rsq.0.3",],
           add.info = TRUE,
           verbose = FALSE)
  {
    ## Read in parameters
    if(is.null(dataset.name))
      dataset.name = basename(dir.name)
    {
      main.signature <- parameter.df[1,"main.signature"]
      correlated.signature <- parameter.df[1,"correlated.signature"]
      sample.number <- parameter.df[1,"sample.number"]
      name.prefix <- parameter.df[1,"name.prefix"]
      main.mean.log <- parameter.df[1,"main.mean.log"]
      main.stdev.log <- parameter.df[1,"main.stdev.log"]
      correlated.stdev.log <- parameter.df[1,"correlated.stdev.log"]
      slope.linear <- parameter.df[1,"slope.linear"]
      main.signature.lower.thres <- parameter.df[1,"main.signature.lower.thres"]
      correlated.signature.lower.thres <- parameter.df[1,"correlated.signature.lower.thres"]
      pearson.r.2.lower.thres <- parameter.df[1,"pearson.r.2.lower.thres"]
      pearson.r.2.higher.thres <- parameter.df[1,"pearson.r.2.higher.thres"]
      min.main.to.correlated.ratio.linear <- parameter.df[1,"min.main.to.correlated.ratio.linear"]
      max.main.to.correlated.ratio.linear <- parameter.df[1,"max.main.to.correlated.ratio.linear"]
      }

    ## Record all the parameters into dataset$parameter
    dataset <- list() # This will contain the data set and the parameters used to generate it
    dataset$parameter <- parameter.df

    ## Set Random Number Generator(RNG) explicitly,
    ## to be compatible with R versions < 3.6.0
    suppressWarnings(RNGkind(sample.kind = "Rounding"))

    ## Set seed
    set.seed(seed)
    seedInUse <- .Random.seed ## Save the seed used so that we can restore the pseudorandom series
    RNGInUse <- RNGkind() ## Save the random number generator (RNG) used


    ## make a directory to store the dataset,
    ## and set the working directory to it.
    if (verbose) {
      cat(paste("Specifying dataset.name as: ", dataset.name,"...\n",sep = ""))
    }
    MustCreateDir(dir.name,overwrite)
    if (verbose) {
      cat(paste("Output folder for this dataset is: ", dir.name,"\n",sep = ""))
    }

    #### Generate exposure matrices for main.signature
    if (verbose) cat("Generating ground-truth exposures...\n")
    dataset$exposure <- GenSBS1SBS5Exposure(main.signature,
                                            correlated.signature,
                                            sample.number,
                                            name.prefix,
                                            main.mean.log,
                                            main.stdev.log,
                                            correlated.stdev.log,
                                            slope.linear,
                                            main.signature.lower.thres,
                                            correlated.signature.lower.thres,
                                            pearson.r.2.lower.thres,
                                            pearson.r.2.higher.thres,
                                            min.main.to.correlated.ratio.linear,
                                            max.main.to.correlated.ratio.linear)

    #### Plot out the scatter plot for the two correlated exposures
    if(add.info){
    if (verbose) cat("Plotting correlation scatterplot for exposures...\n")
    PlotCorrelationScatterplotForExposures(pdf.filename = paste(dir.name,"/scatterplot.pdf",sep = ""),
                                           main.signature = main.signature,
                                           correlated.signature = correlated.signature,
                                           slope.linear,
                                           exposure.counts = dataset$exposure,
                                           xlim = c(0,4),
                                           ylim = c(0,4))
    }

    #### Using the exposure count generate synthetic spectra catalog.
    dataset$spectra <-
      CreateSynCatalogs(sp.sigs[,c(main.signature,correlated.signature)],
                        dataset$exposure)

    if (verbose) cat("Spectra generated.\n")

    #### Output Duke-NUS formatted mutational spectra and exposure.counts
    WriteCatalog(dataset$spectra$ground.truth.catalog,
                 paste0(dir.name,"/ground.truth.syn.catalog.csv"))
    WriteExposure(dataset$exposure,
                  paste0(dir.name,"/ground.truth.syn.exposures.csv"))

    #### Copy ground-truth signatures
    WriteCatalog(dataset$spectra$ground.truth.signatures,
                 paste0(dir.name,"/ground.truth.syn.sigs.csv"))

    #### In debug mode,
    #### Output parameters used for better reproducibility.
    #### In usual cases, we don't output this information,
    #### for passing diff check (NewDiff4SynDataSets)
    if(add.info){
      write.table(t(data.frame(dataset$parameter)),
                  paste0(dir.name,"/parameters.txt"),
                  col.names = F,quote = F,sep = ":\t")
    }

    #### In debug mode,
    #### Save seeds and session information
    #### for better reproducibility.
    #### In usual cases, we don't output this information,
    #### for passing diff check (NewDiff4SynDataSets)
    if(add.info){
      capture.output(sessionInfo(), file = paste0(dir.name,"/sessionInfo.txt")) ## Save session info
      write(x = seedInUse, file = paste0(dir.name,"/seedInUse.txt")) ## Save seed in use to a text file
      write(x = RNGInUse, file = paste0(dir.name,"/RNGInUse.txt")) ## Save seed in use to a text file
    }
    if (verbose) {
      # cat(paste("All result files have been stored in folder ",
      #          dir.name,"\n",sep = ""))
      cat("\n\nData generation finished.\n\n")
    }

  }

#' Function to generate 20 SBS1-SBS5-correlated Synthetic datasets used in testing.
#'
#' This function is a wrapper around \code{\link{CreateSBS1SBS5CorrelatedSyntheticDataOneDataset}}.
#' It will use the default parameters to repeat the results.
#'
#' This function will generate 20 datasets, each with files listed below:
#'
#' ground.truth.syn.catalog.csv: Generated tumor spectra in
#' ICAMS SBS96 CSV format.
#'
#' ground.truth.syn.exposures.csv: Mutation burdens of SBS1 and
#' SBS5 in generated tumor spectra in ICAMS CSV format.
#'
#' ground.truth.syn.sigs.csv: Ground-truth SBS1 and SBS5
#' signatures in ICAMS SBS96 CSV format.
#'
#' parameters.txt: Parameters used to generate the exposures
#' and tumor spectra.
#'
#' scatterplot.pdf: scatterplot illustrating correlation of
#' exposures of two signatures in generated spectra
#'
#' seedInUse.txt, RNGInUse.txt: seed and Random Number Generator
#' used in generation. (For better reproducibility)
#'
#' sessionInfo.txt: information related to R versions, platforms,
#' loaded or imported packages, etc. (For better reproducibility)
#'
#' @param top.level.dir Top-level-folder to place 20 spectra
#' datasets generated by this function.
#' Default: ./ (Current working directory)
#'
#' @param regress.dir If not \code{NULL}, compare the result to
#' the contents of this directory with a \code{diff}.
#'
#' @param overwrite Whether to overwrite
#' (Default: FALSE)
#'
#' @param add.info Whether to generate additional information.
#'
#' You should set it to \code{FALSE} when you want to make a \code{diff}
#' (i.e. \code{regressdir} is not \code{NULL}). This is because
#' Additional information may differ on different OS or R sessions,
#' thus may prevent the dataset from passing the \code{NewDiff4SynDatasets} check.
#' (Default: TRUE)
#'
#' @param unlink Whether to delete temporary dataset folder \code{top.level.dir}.
#' (Set to TRUE for testing)
#'
#' @importFrom ICAMS WriteCatalog
#' @importFrom utils capture.output sessionInfo
#'
#' @export
#'
CreateSBS1SBS5CorrelatedSyntheticData <-
  function(top.level.dir = "./",
           regress.dir = NULL,
           overwrite = FALSE,
           add.info = TRUE,
           unlink = FALSE){

    datasetNames <- rownames(SynSigGen::SBS1SBS5parameter)

    for(datasetName in datasetNames){
      # dataset.name <- datasetName
      dir.name <- paste0(top.level.dir,"/",datasetName,"/")

      ## Call the spectra generation function.
      ## Assign most parameters in a named vector format.
      CreateSBS1SBS5CorrelatedSyntheticDataOneDataset(
        dir.name = dir.name,
        dataset.name = datasetName,
        overwrite = overwrite,
        seed = 1,
        parameter.df = SynSigGen::SBS1SBS5parameter[datasetName,],
        add.info = add.info)
    }

    ## In default mode, the return value will be TRUE
    ## only if the contents in two directories are the same.
    if (!is.null(regress.dir)) {
      diff.result <-
        NewDiff4SynDataSets(newdir         = top.level.dir,
                            regressdirname = regress.dir,
                            unlink         = unlink,
                            verbose        = FALSE,
                            long.diff      = FALSE)
      return(diff.result[1] == "ok")
    }
  }
