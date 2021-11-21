
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SynSigGen

<!-- badges: start -->

[![R-CMD-check](https://github.com/steverozen/SynSigGen/workflows/R-CMD-check/badge.svg?branch=mSigHdp_plus_legacy)](https://github.com/steverozen/SynSigGen/actions?query=workflow%3AR-CMD-check+branch%3AmSigHdp_plus_legacy)

[![R-CMD-check-bioc](https://github.com/steverozen/SynSigGen/workflows/R-CMD-check-bioc/badge.svg?branch=mSigHdp_plus_legacy)](https://github.com/steverozen/SynSigGen/actions?query=workflow%3AR-CMD-check-bioc+branch%3AmSigHdp_plus_legacy)

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
install.packages("remotes")
remotes::install_github(repo = "steverozen/SynSigGen", ref = "master")
```

## Example usage

### 1. Generate synthetic tumor spectra used in paper *The repertoire of mutational signatures in human cancer*

Use functions below to generate 11 spectra datasets used in paper *The
repertoire of mutational signatures in human cancer*
(<https://doi.org/10.1038/s41586-020-1943-3>), published in *Nature*:

``` r
## Users should specify regress.dir = NULL unless while debugging.
PancAdenoCA1000()
ManyTypes2700()
RCCOvary1000()
Create.3.5.40.Abstract()
BladderSkin1000()
Create.2.7a.7b.Abstract()
CreateRandomSyn()
```

The description of 11 data sets are available at section *“Description
of each suite of synthetic data sets”* in Supplementary Note 2 of the
paper.

### Generate synthetic spectra used in paper *“Accuracy of Mutational Signature Software on Correlated Signatures”*

To generate 20 spectra data sets with mutation load of and correlation
between SBS1 and SBS5 varied, use function

``` r
CreateSBS1SBS5CorrelatedSyntheticData()
```

## Notes for functions to generate legacy data sets (1 & 2)

-   The wrapper functions used to generate data sets in [*Nature*
    paper](https://doi.org/10.1038/s41586-020-1943-3) are in R files
    with suffix “\_Nature”.

-   The wrapper functions used to generate data sets for paper on 20
    correlated data sets are in R files with suffix “SBS1SBS5”.

-   These wrapper functions are primarily used to generate legacy data
    sets, as they don’t round the exposures to integers.

-   By contrast, `GenerateSyntheticExposures()` now rounds the exposures
    by default from version 1.0.10.

## Reference manual

<https://github.com/steverozen/SynSigGen/blob/master/data-raw/SynSigGen_1.0.13.pdf>
