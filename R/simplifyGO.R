#' Remove semantically similar GO terms.
#' Based on code by Jamie Soul:
#' https://github.com/soulj/SkeletalVis-Pipeline/blob/master/src/SYBIL%20Systems%20Biology/GOEnrichment.R
#'
#'
#' @param GOID List of GOID's to simplify
#' @param semData output from GOSemSim::godata
#' @importFrom GO.db GOBPCHILDREN
#' @importFrom GOSemSim goSim
#' @importFrom dplyr filter
#' @return vector of reduced GO terms

simplifyGO <- function(GOID, semData) {
  semSim <- expand.grid(unique(GOID),
                        unique(GOID), stringsAsFactors = F) %>%
    dplyr::filter(.data$Var1 != .data$Var2)
  semSim$sim <-
    apply(semSim, 1, function(i)
      GOSemSim::goSim(i[1], i[2], semData = semData, measure = "Wang"))

  semSim[is.na(semSim)] <- 0
  semSim <- semSim[!is.na(semSim$sim), ]
  semSim <- semSim[order(semSim$sim , decreasing = T), ]

  #get the similar term pairs
  semSim <- semSim[semSim$Var1 != semSim$Var2, ]
  semSim <- semSim[semSim$sim > 0.4, ]


  #mark high fequency terms
  semSim$remove <- apply(semSim, 1, function(x) {
    if (x[1] %in% highFreqTerms) {
      return(x[1])
    }
    if (x[2] %in% highFreqTerms) {
      return(x[2])
    } else {
      return(NA)
    }
  })
  remove <- na.omit(semSim$remove)
  semSim <- semSim[is.na(semSim$remove), ]

  if (nrow(semSim) == 0) {
    NA
  } else{
    for (i in 1:nrow(semSim)) {
      Var1 <- semSim[i, "Var1"]
      Var2 <- semSim[i, "Var2"]
      if (Var2 %in% as.list(GO.db::GOBPCHILDREN)[[Var1]]) {
        remove <- c(remove, Var2)
        next
      } else if (Var1 %in% as.list(GO.db::GOBPCHILDREN)[[Var2]])
        remove <- c(remove, Var1)
      next
    }
  }

  flt <- unique(GOID[!(GOID %in% remove)])

  return(flt)
}
