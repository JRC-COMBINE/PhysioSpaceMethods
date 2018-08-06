#' @title Mapping one sample into a physiological-space
#'
#' @description SingleThreadOfPhysioCalc is an internal function of calculatePhysioMap,
#' computing the main mapping. We don't recommend the use of
#' SingleThreadOfPhysioCalc outside of calculatePhysioMap.
#'
#' @import progress parallel
#'
#' @param INPT A sample (column) out of InputData.
#' @param Space The space in which the 'InputData' will be mapped. Just as 'InputData',
#' it should be a matrix with genes as rows and samples as columns, with corresponding Entrez Gene
#' IDs in 'rownames' of the matrix, and name of each axis of the space written in 'colnames'.
#' Rows of InputData and Space must match.
#' @param GenesRatio The ratio of gene expression values to be considered in the calculation. In high dimensional omics data,
#' signal to noise ratio has a direct relation with the relative magnitude of expressions. We aim to remove the noisy genes,
#' hence we only keep the "GenesRatio*100" percent highest and lowest gene expression values of each sample.
#' GenesRatio should be a numerical value between 0 and 1. Default value is 0.05.
#' @param TTEST Logical value indicating if t.test should be done in place of the default wilcoxon rank-sum test (more info can be found
#' in the original PhysioSpace: Lenz et. al., PLOS One 2013). Using t.test will speed up calculations. Default value is FALSE.
#' @param STATICResponse Logical value indicating if 'statistic' should be returned rather than the default 'signed p value'.
#' Default value is FALSE.
#' @param NGenes Number of genes (rows) in Space.
#' @param pb Progress bar made by progress::progress_bar$new.
#'
#' @return Vector of mapped 'INPT' values in 'Space'. Mapped values are signed p value when STATICResponse==FALSE,
#' and are 'statistic' value when STATICResponse==TRUE (more info can be found in the original PhysioSpace:
#' Lenz et. al., PLOS One 2013).
#'
# #' @export #Not exporting since it's an internal function
singleThreadOfPhysioCalc <- function(INPT, Space, GenesRatio, NGenes, STATICResponse, pb, TTEST){
  if (!is.null(GenesRatio)) {
    numGenes = round(NGenes * GenesRatio)
  } else{
    stop("Don't have Number of genes to compare??!?!")
  }
  ordDiff = order(INPT)
  if (!is.null(numGenes)) {
    iplus = ordDiff[seq.int(from = (NGenes - numGenes + 1), to = NGenes)]
    iminus = ordDiff[seq_len(numGenes)]
  }
  pb$tick()
  apply(X = Space, MARGIN = 2, FUN = if(TTEST) tTestWrapper else wilTestWrapper,
        iplus=iplus, iminus = iminus, STATICResponse = STATICResponse)
}
