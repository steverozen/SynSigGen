# Warning, some signatures present in only one sample, dropping: SBS43 SBS53
SBS96.real.exposures <- PCAWG7::exposure$PCAWG$SBS96
SBS96.sig.params <-
  GetSynSigParamsFromExposures(exposures = SBS96.real.exposures,
                               verbose = 1,
                               distribution = "neg.binom")

DBS78.real.exposures <- PCAWG7::exposure$PCAWG$DBS78
DBS78.sig.params <-
  GetSynSigParamsFromExposures(exposures = DBS78.real.exposures,
                               verbose = 1,
                               distribution = "neg.binom")

ID.real.exposures <- PCAWG7::exposure$PCAWG$ID
ID.sig.params <-
  GetSynSigParamsFromExposures(exposures = ID.real.exposures,
                               verbose = 1,
                               distribution = "neg.binom")

sig.params <- list(SBS96 = SBS96.sig.params,
                   DBS78 = DBS78.sig.params,
                   ID = ID.sig.params)
usethis::use_data(sig.params)
