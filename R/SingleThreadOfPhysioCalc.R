#' @title Mapping one sample into a physiological-space
#'
#' @description SingleThreadOfPhysioCalc is an internal function of calculatePhysioMap,
#' computing the main mapping. It's not recommended to use SingleThreadOfPhysioCalc by the user outside of
#' calculatePhysioMap.
#'
#' @import progress parallel
#'
# #' @export #Not exporting since it's an internal function
singleThreadOfPhysioCalc <- function(INPT, Space, NSamples, GenesRatio, NGenes, STATICResponse, pb, TTEST){
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
