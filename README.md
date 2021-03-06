
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SynSigGen

<!-- badges: start -->

[![R-CMD-check](https://github.com/steverozen/SynSigGen/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/steverozen/SynSigGen/actions?query=workflow%3AR-CMD-check+branch%3Amaster)

[![R-CMD-check-bioc](https://github.com/steverozen/SynSigGen/workflows/R-CMD-check-bioc/badge.svg?branch=master)](https://github.com/steverozen/SynSigGen/actions?query=workflow%3AR-CMD-check-bioc+branch%3Amaster)

<!-- badges: end -->

Synthetic (Mutational) Signature Generation (**`SynSigGen`**)

## Purpose

Create catalogs of synthetic mutational spectra for assessing the
performance of mutational signature analysis programs. ‘SynSigGen’
stands for Generation of Synthetic Signatures and Spectra.

## Installation

Before installation, prerequisites in
[Bioconductor](https://www.bioconductor.org/) needs to be installed:

``` r
install.packages("BiocManager")
BiocManager::install(c("Biostrings","BSgenome","GenomeInfoDb","GenomicRanges"))
```

Install from GitHub with the R command line:

``` r
install.packages("devtools")
devtools::install_github("steverozen/SynSigGen")
```

## Example usage

### Generate synthetic tumor spectra used in paper *The repertoire of mutational signatures in human cancer*

Use functions below to generate 11 spectra datasets used in paper *The
repertoire of mutational signatures in human cancer*
(<https://doi.org/10.1038/s41586-020-1943-3>):

``` r
PancAdenoCA1000()
ManyTypes2700()
RCCOvary1000()
Create.3.4.40.Abstract()
BladderSkin1000()
Create.2.7a.7b.Abstract()
CreateRandomSyn()
```

### Generate synthetic spectra used in paper *“Performance of Mutational Signature Software on Correlated Signatures”*

To generate 20 spectra data sets with mutation load of and correlation
between SBS1 and SBS5 varied, use function

``` r
CreateSBS1SBS5CorrelatedSyntheticData()
```

## Reference manual

<https://github.com/steverozen/SynSigGen/blob/master/data-raw/SynSigGen_1.0.6.pdf>
