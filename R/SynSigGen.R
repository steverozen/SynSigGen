#' @title SynSigGen
#'
#' @description Create catalogs of synthetic mutational spectra for
#' assessing the performance of mutational-signature analysis programs.
#'
#' @section Overview:
#'
#' The main focus is generating synthetic catalogs of mutational
#' spectra (mutations in tumors) based on known mutational signature
#' profiles and software-inferred exposures (software's estimate on 
#' number of mutations induced by mutational signatures in tumors) 
#' in the PCAWG7 data. We call this kind of synthetic data broadly 
#' "reality-based" synthetic data. 
#' The package also has a set of functions that generate
#' random mutational signature profiles and then create synthetic 
#' mutational spectra based on these random signature profiles. We
#' call this kind of synthetic data "random" synthetic data, while
#' pointing out that much depends on the distributions from which
#' the random signature profiles and attributions are generated.
#'
#' @section Workflow for generating "reality-based" synthetic mutational spectra:
#'
#' Typical workflow for generating catalogs of "reality-based" synthetic
#' mutational spectra is as follows. \enumerate{
#'
#' \item Input (based on SignatureAnalyzer or SigProfiler analysis of PCAWG tumors)
#'   \code{E}, matrix of software-inferred exposures of mutational signatures (signatures x samples)
#'   \code{S}, mutational signature profiles (mutation types x signatures)
#'
#' \item Obtain distribution parameters from software-inferred exposures \preformatted{
#'   P <- GetSynSigParamsFromExposures(E, ...)
#' }
#'
#' \item Generate exposures for synthetic mutational spectra based on \code{P} \preformatted{
#'   synthetic.exposures <- GenerateSyntheticExposures(P, ...)
#' }
#'
#' \item Generate synthetic mutational spectra by multiplying \code{S} and \code{synthetic.exposures},
#' and round the product to the nearest unit: \preformatted{
#'   synthetic.spectra <- CreateAndWriteCatalog(S, synthetic.exposures, ...)
#' }
#'
#' }
#'
#' @section Workflow for generating "random" synthetic mutational spectra:
#'
#' The top-level function for generating "random" synthetic mutational spectra is 
#' \code{\link{CreateRandomSyn}}. It adopts the following steps to generate
#' catalogs of "random" synthetic mutational spectra. \enumerate{
#'
#' \item Create random mutational signature profiles: \preformatted{
#'   S <- CreateRandomMutSigProfiles(...)
#' }
#'
#' \item Generate distribution parameters for exposures of random signatures: \preformatted{
#'   P <- CreateMeanAndStdevForSigs(sig.names = colnames(S),...)
#' }
#'
#' \item Create exposures for mutational signatures based on \code{P} and other 
#' parameters: \preformatted{
#'   synthetic.exposures <- CreateRandomExposures(sigs = S, per.sig.mean.and.sd = P)
#' }
#'
#' \item Generate synthetic mutational spectra by multiplying \code{S} and \code{synthetic.exposures}
#' and round the product to the nearest unit: \preformatted{
#'   synthetic.spectra <- NewCreateAndWriteCatalog(S, synthetic.exposures, ...)
#' }
#'
#' }
#'
#' @section Workflow for generating "SBS1-SBS5-correlated" synthetic mutational spectra
#'
#' Typical workflow for generating catalogs of "SBS1-SBS5-correlated" synthetic
#' mutational spectra is as follows. \enumerate{
#'
#'
#'
#' }
#'
#'
#'
#'
#' @docType package
#' @name SynSigGen
#'
NULL
