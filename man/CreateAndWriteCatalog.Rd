% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CreateSynData.R
\name{CreateAndWriteCatalog}
\alias{CreateAndWriteCatalog}
\title{Create and write a mutational spectra catalog}
\usage{
CreateAndWriteCatalog(
  sigs,
  exp,
  dir = NULL,
  write.cat.fn = ICAMS::WriteCatalog,
  extra.file.suffix = "",
  overwrite = FALSE,
  my.dir = NULL
)
}
\arguments{
\item{sigs}{Signatures to use.}

\item{exp}{(Synthetic) exposures.}

\item{dir}{Deprecated, maintained only to avoid
breaking old code. A subdirectory based on
the deprecated global variable \code{\link{OutDir}}.}

\item{write.cat.fn}{Function to write catalogs \strong{or}
spectra to files.}

\item{extra.file.suffix}{Extra string to put before ".csv".}

\item{overwrite}{If TRUE, overwrite existing directory; useful for
debugging / testing.}

\item{my.dir}{The directory in which to write the catalog
and several additional files.}
}
\value{
Invisibly, the generated catalog.
}
\description{
Create and write a mutational spectra catalog
}
\details{
Create a file with the catalog \code{syn.data.csv}
 and writes \code{sigs} to \code{input.sigs.csv}.
}
