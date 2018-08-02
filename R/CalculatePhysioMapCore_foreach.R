#' @title Mapping new data into a physiological-space using foreach
#'
#' @description calculatePhysioMapCore_foreach is an internal function of calculatePhysioMap,
#' computing the main mapping. It's not recommended to use calculatePhysioMapCore_foreach by the user outside of
#' calculatePhysioMap.
#'
#' @import progress doParallel foreach parallel
#'
# #' @export #Not exporting since it's an internal function
calculatePhysioMapCore_foreach <- function(InputData, Space, NSamples, GenesRatio, NGenes, STATICResponse, pb, TTEST){
  foreach(SAMPEL=seq_len(NSamples), .combine='cbind', .final=as.matrix) %dopar% {
    tempDiff <- InputData[, SAMPEL]
    if (!is.null(GenesRatio)) {
      numGenes = round(NGenes * GenesRatio)
    } else{
      stop("Don't have Number of genes to compare??!?!")
    }
    ordDiff = order(tempDiff)
    if (!is.null(numGenes)) {
      iplus = ordDiff[seq.int(from = (NGenes - numGenes + 1), to = NGenes)]
      iminus = ordDiff[seq_len(numGenes)]
    }
    pb$tick()
    apply(X = Space, MARGIN = 2, FUN = if(TTEST) tTestWrapper else wilTestWrapper,
          iplus=iplus, iminus = iminus, STATICResponse = STATICResponse) # this apply can be written as simple for loop also, results
    ## seems to be the same, but writing a nested foreach loop may speed up stuff -> will have to try later
  }
}
