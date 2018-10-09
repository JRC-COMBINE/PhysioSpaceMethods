#' @title Mapping one sample into a physiological-space
#'
#' @description .singleThreadOfPhysioCalc is an internal function
#' of calculatePhysioMap, computing the main mapping. We don't
#' recommend the use of .singleThreadOfPhysioCalc outside of
#' calculatePhysioMap.
#'
#' @import progress parallel
#'
#'
#' @inheritParams calculatePhysioMap
#'
#' @param INPT A sample (column) out of InputData.
#' @param NGenes Number of genes (rows) in Space.
#' @param pb Progress bar made by progress::progress_bar$new.
#'
#' @return Vector of mapped 'INPT' values in 'Space'.
#' Mapped values are signed p value when STATICResponse==FALSE,
#' and are 'statistic' value when STATICResponse==TRUE
#' (more info can be found in the original PhysioSpace:
#' Lenz et. al., PLOS One 2013).
#'
# #' @export #Not exporting since it's an internal function
.singleThreadOfPhysioCalc <- function(INPT, Space, GenesRatio,
                                        NGenes, STATICResponse, pb, TTEST){
    if (!is.null(GenesRatio)) {
        numGenes = round(NGenes * GenesRatio)
    } else{
        stop("'GenesRatio' must be provided.")
    }
    ordDiff = order(INPT)
    if (!is.null(numGenes)) {
        iplus = ordDiff[seq.int(from = (NGenes - numGenes + 1), to = NGenes)]
        iminus = ordDiff[seq_len(numGenes)]
    }
    pb$tick()
    apply(X = Space, MARGIN = 2,
            FUN = if(TTEST) .tTest else .wilTest,
            iplus=iplus, iminus = iminus, STATICResponse = STATICResponse)
}
