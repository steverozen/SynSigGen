## SynSigGen 1.2.1
### Fixed
* Changed functions that were moved from mSigAct to mSigTools.
* Minor documentation fixes and fixed use of == to test class.

### Added
* Added new internal function `GetNonZeroNoisySample` to add noise to one sample until the final spectrum has mutation > 0.

### Updated
* Updated exported function `AddNoise` to make sure the samples in the final noisy spectra all have mutations > 0.

<br>

## SynSigGen 1.2.0
### Fixed
* Fixed a bug in `GenerateSyntheticExposures` when generating only one synthetic tumor
with only one signature.

<br> 

## SynSigGen 1.1.1
### Fixed
* Fixed a bug in function `GenerateNoisyTumors` to only write signatures that
are present in exposures to CSV file.

### Added
* Added new dependency packages `cosmicsig` and `mSigAct`.

<br> 

## SynSigGen 1.1.0
### Backward compatibility fix

Due to fix in function `GenerateSynExposureOneSample` to round the mutations
due to each signature in version 1.0.9, the wrapper functions to generate 
data sets for Nature paper and SBS1-SBS5 paper fail to reproduce the data sets
used by these papers.

In SynSigGen 1.1.0, we disabled the rounding of mutations due to each signature
when running generator functions for legacy data sets, so that these legacy data
sets can be reproduced.

### Updated
* Updated documentation for argument `sig.params` in function `GenerateSyntheticExposures`.

### Added
* Added two fields `URL` and `BugReports` in DESCRIPTION.

<br> 

## SynSigGen 1.0.13
### Removed
* Removed exported data `ID.MMR.params` because the list of MSI-H tumors is not
complete when generating the old data.

<br> 

## SynSigGen 1.0.12
### Update 
* Updated exported function `MergeExposures` to sort the signature ids of the 
merged exposure.

* Suppressed warning messages that came from internal function `NumFromId`.

### Bug fix
* Fixed bug in internal function `GetMutationType` to check all the signature
names to determine the mutation type.

<br> 

## SynSigGen 1.0.11
### **Important** Update 
* Updated exported function `GetSynSigParamsFromExposures` and internal function
`SynSigParamsOneSignature` to only use the `size` parameter from empirical
signature parameters if there is only one data to fit the negative binomial
distribution. The `mu` parameter will be using the original mutation count.

<br> 

## SynSigGen 1.0.10
### Updated 
* Updated function `GenerateSyntheticTumors` to show informative messages
about regenerating parameters from the synthetic exposures and compare with that
from real exposures if `verbose > 0`.

### Added
* Added new argument `tumor.marker.name` in function `GenerateSyntheticTumors`.

* Added new argument `sig.params` in exported functions
`GetSynSigParamsFromExposures` and `GenerateSyntheticTumors`.

* Added two new exported functions `GenerateListOfSigParams` and
`GenerateSyntheticTumorsFromSigParams`.

<br> 

## SynSigGen 1.0.9
### Updated
* Updated function  `GenerateSynExposureOneSample` to round the synthetic exposure to make
it biologically reasonable.

* Updated function `GetSynSigParamsFromExposures` not to drop rare signatures
that are present in only one sample when `distribution` is `neg.binom`. Instead,
this function will use the empirical signature parameters from all cancer types.

* Updated function `GenerateSyntheticTumors` to enable generating different number 
of synthetic tumors for each cancer type.

* Updated function `GetSynSigParamsFromExposures` to drop signatures which do not
have empirical signature parameters from all cancer types.

* Updated function `GenerateSyntheticTumors` and `GenerateNoisyTumors` to append
mutation type information to the names of files generated.

* Updated function `GenerateSyntheticExposures` to remove signatures that have zero
exposure in the synthetic data.

### Added
* Added default value 100 for argument `n.binom.size` in function `AddNoise`. 

* Added new exported data `signature.params`.

* Added new exported data `ID.MMR.params`.

* Added two new internal functions `NumFromId` and `SortSigId`.

### Fixed
* Fixed a bug in function `GenerateSynExposureOneSample` to round the mutations
due to each signature.

<br> 

## SynSigGen 1.0.8
### Added 
* Added new exported function  `GenerateNoisyTumors`.

* Minor change to function `GenerateSyntheticTumors` for always writing the
signature names in the parameters CSV file in case some signature is not present
in the synthetic data.

<br> 

## SynSigGen 1.0.7
### Added 
* Added new argument `distribution` in functions `SynSigParamsOneSignature` and
`GetSynSigParamsFromExposures` to enable fitting the exposures using negative
binomial distribution.

* Added new argument `distribution` and `sig.params` in function
`GenerateSynExposureOneSample` to enable generating synthetic exposures using
negative binomial distribution.

* Added new argument `distribution` in function `GenerateSyntheticExposures` to
enable generating synthetic exposures using negative binomial distribution.

* Added new exported function `GenerateSyntheticTumors` which enables generating 
synthetic tumors using negative binomial distribution.

* Added new argument `verbose` in function `GenerateSyntheticTumors` which can
cat various messages if value greater than 0.

* Added new argument `cancer.type` in function `GetSynSigParamsFromExposures`.

### Fixed
* Fixed a bug in function `SynSigParamsOneSignature` to generate parameter
estimates (negative binomial distribution) for one signature only using those
tumors which have mutations for this particular signature.

### Updated
* Updated function `SynSigParamsOneSignature` to return NA for the parameter
estimates (negative binomial distribution) for one signature if there is only
one tumor which has mutation for this particular signature.

<br> 

## SynSigGen 1.0.6
### Added 
* Added data object "SBS1SBS5datasetNames" as names of SBS1-SBS5 spectra datasets generated by CreateSBS1SBS5CorrelatedSyntheticData()

* Added new argument `sig.matrix` in functions `GenerateSynExposureOneSample` and `GenerateSyntheticExposures` so that we can resample the synthetic exposures to make
sure the total mutations for the reconstructed catalog is not zero.

### Fixed 
* Fixed a bug in `GetSynSigParamsFromExposures` when the exposures only have one row that are
not all zeros.

## SynSigGen 1.0.5
### Simplified
* Simplified folder structure of SBS1-SBS5-correlated dataset,
generated by `CreateSBS1SBS5CorrelatedSyntheticData()` - 
The files of synthetic datasets are under the dataset folders, and no more
`sp.sp` folders.

## SynSigGen 1.0.4
### Last version generates SBS1-SBS5-correlated dataset with two-layer directory  structure:

* Name of dataset (e.g. `S.0.1.Rsq.0.1`)
* Folder `sp.sp`, means that two ground-truth signatures (`SBS1` and `SBS5`) are from results of SigProfiler in (https://doi.org/10.1038/s41586-020-1943-3)
* Files of synthetic datasets.

### Clarified
* Clarified package `README`
* Clarified package documentation
* Clarified function documentation

### Fixed typo 
* Renamed function `Create.3.4.40.Abstract()` to `Create.3.5.40.Abstract()`

## SynSigGen 1.0.3
* Removed `R/SigSimilarity.R` (functionality now in `ICAMSxtra`)
  and adjusted code in this package.

## SynSigGen 1.0.2.9001

### Added
* Added function `AddNoise`.
* Added `biocViews` as a field name in package [DESCRIPTION](https://github.com/steverozen/SynSigGen/blob/master/DESCRIPTION).
* Created reference manual.
