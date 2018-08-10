#' @title Mapping new data into a physiological-space using foreach
#'
#' @description calculatePhysioMapCore_foreach is an internal function of
#' calculatePhysioMap, computing the main mapping using foreach form
#' foreach package. We don't recommend the use of
#' calculatePhysioMapCore_foreach outside of calculatePhysioMap.
#'
#' @import progress doParallel foreach parallel
#'
#' @inheritParams calculatePhysioMap
#'
#' @param NSamples Number of samples (columns) in InputData.
#' @param NGenes Number of genes (rows) in InputData or Space.
#' @param pb Progress bar made by progress::progress_bar$new.
#'
#' @return Matrix of mapped 'InputData' values in 'Space', with
#' rows corresponding to axes of 'Space' and columns representing
#' samples in 'InputData'. Mapped values are signed p value when
#' STATICResponse==FALSE, and are 'statistic' value when
#' STATICResponse==TRUE (more info can be found in the original
#' PhysioSpace paper: Lenz et. al., PLOS One 2013).
#'
# #' @export #Not exporting since it's an internal function
calculatePhysioMapCore_foreach <- function(InputData, Space,
                                            NSamples, GenesRatio,
                                            NGenes, STATICResponse,
                                            pb, TTEST){
    foreach(SAMPEL=seq_len(NSamples),
            .combine='cbind', .final=as.matrix) %dopar% {
                tempDiff <- InputData[, SAMPEL]
                if (!is.null(GenesRatio)) {
                    numGenes = round(NGenes * GenesRatio)
                } else{
                    stop("Don't have Number of genes to compare??!?!")
                }
                ordDiff = order(tempDiff)
                if (!is.null(numGenes)) {
                    iplus = ordDiff[seq.int(from = (NGenes - numGenes + 1),
                                            to = NGenes)]
                    iminus = ordDiff[seq_len(numGenes)]
                }
                pb$tick()
                apply(
                    X = Space,
                    MARGIN = 2,
                    FUN = if (TTEST)
                        tTestWrapper
                    else
                        wilTestWrapper,
                    iplus = iplus,
                    iminus = iminus,
                    STATICResponse = STATICResponse
                ) # this apply can be written as simple for loop also, results
                ## seems to be the same, but writing a nested foreach loop may
                # speed up stuff -> will have to try later
            }
}
