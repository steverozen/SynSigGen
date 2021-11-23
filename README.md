
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SynSigGen

<!-- badges: start -->

[![R-CMD-check](https://github.com/steverozen/SynSigGen/workflows/R-CMD-check/badge.svg?branch=1.1.0-branch)](https://github.com/steverozen/SynSigGen/actions?query=workflow%3AR-CMD-check+branch%3A1.1.0-branch)

[![R-CMD-check-bioc](https://github.com/steverozen/SynSigGen/workflows/R-CMD-check-bioc/badge.svg?branch=1.1.0-branch)](https://github.com/steverozen/SynSigGen/actions?query=workflow%3AR-CMD-check-bioc+branch%3A1.1.0-branch)

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
(<https://doi.org/10.1038/s41586-020-1943-3>), published in *Nature*.
The data sets are available at
[Synapse](https://www.synapse.org/#!Synapse:syn18497223):

``` r
# Users should specify regress.dir = NULL unless for comparison
# with original data sets.
# 
# Compare tools (e.g. BeyondCompare, Meld) is recommended 
# over specifying regress.dir,
# because the latter might raise an error even when query 
# and original data sets are identical.
#
# Users should specify top.level.dir to the destination folder
# for data sets. Otherwise default paths will be used.
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

The original data sets are available at
[Zenodo](https://doi.org/10.5281/zenodo.2636980).

## Notes for functions to generate legacy data sets (1 & 2)

-   The wrapper functions used to generate data sets in [*Nature*
    paper](https://doi.org/10.1038/s41586-020-1943-3) are in R files
    with suffix “\_Nat”.

-   The wrapper functions used to generate data sets for paper on 20
    correlated data sets are in file “CreateSynSBS1SBS5Correlated.R”.

-   These wrapper functions are primarily used to generate legacy data
    sets, as they don’t round the exposures to integers.

-   By contrast, `GenerateSyntheticExposures()` now rounds the exposures
    by default from version 1.0.10.

## Reference manual

<https://github.com/steverozen/SynSigGen/blob/mSigHdp_plus_legacy/data-raw/SynSigGen_1.1.0.pdf>
