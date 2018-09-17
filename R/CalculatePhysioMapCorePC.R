#' @title Mapping new data into a physiological-space using parCapply
#' of parallel package
#'
#' @description calculatePhysioMapCorePC (PC: parCapply) is an internal
#' function of calculatePhysioMap, computing the main mapping using
#' pacCapply of parallel package. We don't recommend the use of
#' calculatePhysioMapCorePC outside of calculatePhysioMap.
#'
#' @import progress parallel
#'
#' @inheritParams calculatePhysioMap
#'
#' @param NSamples Number of samples (columns) in InputData.
#' @param NGenes Number of genes (rows) in InputData or Space.
#' @param pb Progress bar made by progress::progress_bar$new.
#' @param cl Cluster made by parallel::makeCluster.
#'
#' @return  Matrix of mapped 'InputData' values in 'Space', with
#' rows corresponding to axes of 'Space' and columns representing
#' samples in 'InputData'. Mapped values are signed p value when
#' STATICResponse==FALSE, and are 'statistic' value when
#' STATICResponse==TRUE (more info can be found in the original
#' PhysioSpace paper: Lenz et. al., PLOS One 2013).
#'
# #' @export #Not exporting since it's an internal function
calculatePhysioMapCorePC <- function(InputData, Space,
                                                NSamples, GenesRatio,
                                                NGenes, STATICResponse,
                                                pb, TTEST, cl){
    matrix(parCapply(cl = cl, x = InputData, FUN = singleThreadOfPhysioCalc,
                        Space=Space, GenesRatio=GenesRatio, NGenes=NGenes,
                        STATICResponse=STATICResponse, pb=pb, TTEST=TTEST),
            ncol = ncol(InputData))
}
