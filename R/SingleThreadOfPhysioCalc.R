#' @title Mapping one sample into a physiological-space
#'
#' @description .singleThreadOfPhysioCalc is an internal function
#' of calculatePhysioMap, computing the main mapping. We don't
#' recommend the use of .singleThreadOfPhysioCalc outside of
#' calculatePhysioMap().
#'
#' @import progress
#'
#' @inheritParams calculatePhysioMap
#'
#' @param InputData A matrix, SummarizedExperiment object or a list,
#' based on the gene expression data (or any other type of high
#' dimensional data, e.g. protein abundance, SNP, Methylation, etc.),
#' to be analysed. InputData has to have a specific format to be
#' properly analysed, these requirements are thoroughly explained
#' in the 'Details' section of calculatePhysioMap() function.
#' @param SampleNum A sample (column) number of InputData.
#' @param NGenes Number of genes (rows) in Space.
#' @param pb Progress bar made by progress::progress_bar$new.
#'
#' @return Vector of mapped 'InputData[,SampleNum]' values in 'Space'.
#' Mapped values are signed p value when STATICResponse==FALSE,
#' and are 'statistic' value when STATICResponse==TRUE
#' (more info can be found in the original PhysioSpace:
#' Lenz et. al., PLOS One 2013).
#'
# #' @export #Not exporting since it's an internal function
.singleThreadOfPhysioCalc <- function(SampleNum, InputData, Space, GenesRatio,
                                        NGenes, STATICResponse, pb, TTEST){
    if (!is.null(GenesRatio)) {
        numGenes = round(NGenes * GenesRatio)
    } else{
        stop("'GenesRatio' must be provided.")
    }
    ordDiff = order(InputData[,SampleNum])
    if (!is.null(numGenes)) {
        iplus = ordDiff[seq.int(from = (NGenes - numGenes + 1), to = NGenes)]
        iminus = ordDiff[seq_len(numGenes)]
    }
    pb$tick()
    apply(X = Space, MARGIN = 2,
            FUN = if(TTEST) .tTest else .wilTest,
            iplus=iplus, iminus = iminus, STATICResponse = STATICResponse)
}
