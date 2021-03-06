% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CreateSynSBS1SB5Correlated.R
\name{GenSBS1SBS5Exposure}
\alias{GenSBS1SBS5Exposure}
\title{Generate correlated exposures for multiple tumors

Wrapper function around \code{\link{GenSBS1SBS5ExposureOneTumor}}:
A function to generate exposure of two correlated signatures
(Example: SBS1 and SBS5) for \code{sample.number} (e.g. 500) synthetic tumors.

NOTE: \code{pearson.r.2.lower.thres} and \code{pearson.r2.higher.thres}
are used to constraint the Pearson's R^2 of mutation burdens of two signatures
in multiple tumors.}
\usage{
GenSBS1SBS5Exposure(
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
  pearson.r.2.higher.thres = 1,
  min.main.to.correlated.ratio.linear = 1/3,
  max.main.to.correlated.ratio.linear = Inf
)
}
\arguments{
\item{main.signature}{Name of a signature with smaller variance
in the log10 space. (Default: "SBS5")}

\item{correlated.signature}{Name of a signature with larger variance
in the log10 space. (Default: "SBS1")}

\item{sample.number}{Number of tumors whose mutation burdens
will be generated. (Default: 500)}

\item{name.prefix}{Prefix of tumor name.
(Default: \code{"TwoCorreSigsGen"})
By default, the name of tumors to be created will be:
TwoCorreSigGen::1, TwoCorreSigGen::2, TwoCorreSigGen::3...}

\item{main.mean.log}{Mean of log10(mutation burden of \code{main.signature})}

\item{main.stdev.log}{Standard deviation of log10(mutation burden
of \code{main.signature})}

\item{correlated.stdev.log}{Contribute to part of the standard deviation of
log10(mutation burden of correlated.signature). In this script, the s.d. of
log10(mutation burden of correlated.signature)
= main.stdev.log + correlated.stdev.log}

\item{slope.linear}{Average ratio of mutation burden of \code{correlated.signature}
over mutation burden of \code{main.signature}}

\item{main.signature.lower.thres}{Minimum mutation burden
(number of mutations) induced by \code{main.signature} in each tumor.}

\item{correlated.signature.lower.thres}{Minimum mutation burden
(number of mutations) induced by \code{correlated.signature} in each tumor.}

\item{pearson.r.2.lower.thres}{Minimum Pearson's R^2 of mutation burdens
of two signatures in \code{sample.number} tumors.}

\item{pearson.r.2.higher.thres}{Maximum Pearson's R^2 of mutation burdens
of two signatures in \code{sample.number} tumors.}

\item{min.main.to.correlated.ratio.linear}{Minimum ratio of
\code{main.signature} over mutation burden of
\code{correlated.signature} in each tumor.}

\item{max.main.to.correlated.ratio.linear}{Maximum ratio of
\code{main.signature} over mutation burden of
\code{correlated.signature} in each tumor.}
}
\description{
Generate correlated exposures for multiple tumors

Wrapper function around \code{\link{GenSBS1SBS5ExposureOneTumor}}:
A function to generate exposure of two correlated signatures
(Example: SBS1 and SBS5) for \code{sample.number} (e.g. 500) synthetic tumors.

NOTE: \code{pearson.r.2.lower.thres} and \code{pearson.r2.higher.thres}
are used to constraint the Pearson's R^2 of mutation burdens of two signatures
in multiple tumors.
}
