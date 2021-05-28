MMR.samples <- c("Kidney-RCC::SP102897", "Skin-Melanoma::SP124384",
                 "Liver-HCC::SP107012", "Liver-HCC::SP98845",
                 "Ovary-AdenoCA::SP102133")
ID.exposures <- PCAWG7::exposure$PCAWG$ID
MMR.sample.indices.exposure <- sapply(MMR.samples, FUN = function(x) {
  sample.name <- x
  index <- grep(pattern = sample.name, x = colnames(ID.exposures))
})
ID.exposures.MMR <- ID.exposures[, MMR.sample.indices.exposure, drop = FALSE]

ID.MMR.params <-
  GetSynSigParamsFromExposuresOld(exposures = ID.exposures.MMR,
                                  verbose = 1,
                                  distribution = "neg.binom")

usethis::use_data(ID.MMR.params)
