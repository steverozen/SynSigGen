#' Create a specific synthetic data set of 2,700 tumors.
#'
#' @keywords internal

Create.syn.many.types <-
  function(regress = FALSE, seed = NULL, unlink = FALSE) {
    stop("This function is no longer used")

    if (is.null(seed)) {
      suppressWarnings(RNGkind(sample.kind = "Rounding"))
      # For compatibility with R < 3.6.0
      set.seed(191906)
      top.level.dir <- "tmp.syn.many.types"

    } else {
      set.seed(seed)
      top.level.dir <- paste0("../2700.tumors.seed.", seed)
    }
    num.syn.tumors <- 300 # number of tumor of each type

    cancer.types <- c("Bladder-TCC",    "Eso-AdenoCA",
                      "Breast-AdenoCA", "Lung-SCC",
                      "Kidney-RCC",     "Ovary-AdenoCA",
                      "Bone-Osteosarc", "Cervix-AdenoCA",
                      "Stomach-AdenoCA")
    retval <-
      CreateMixedTumorTypeSyntheticData(
        top.level.dir = top.level.dir,
        cancer.type.strings = cancer.types,
        num.syn.tumors = num.syn.tumors,
        overwrite = TRUE
      )

    if (regress) {
      diff.result <- Diff4SynDataSets("syn.many.types", unlink = unlink)
      if (diff.result[1] != "ok") {
        message("\nThere was a difference, investigate\n",
                paste0(diff.result, "\n"))
      } else {
        message("\nok\n")
      }
    }

    invisible(retval)
  }

for.ludmil.2019.08.16 <- function() {
  seeds <- sample(10000, size = 9)
  for (seed in seeds) {
    Create.syn.many.types(seed = seed)
  }
}

