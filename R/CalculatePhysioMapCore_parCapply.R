#' @title Mapping new data into a physiological-space using parCapply of parallel package
#'
#' @description calculatePhysioMapCore_parCapply is an internal function of calculatePhysioMap,
#' computing the main mapping. It's not recommended to use calculatePhysioMapCore_parCapply by the user outside of
#' calculatePhysioMap.
#'
#' @import progress parallel
#'
# #' @export #Not exporting since it's an internal function
calculatePhysioMapCore_parCapply <- function(InputData, Space, NSamples, GenesRatio, NGenes, STATICResponse, pb, TTEST, cl){
  matrix(parCapply(cl = cl, x = InputData, FUN = singleThreadOfPhysioCalc, Space=Space,
                   NSamples=NSamples, GenesRatio=GenesRatio, NGenes=NGenes,
                   STATICResponse=STATICResponse, pb=pb, TTEST=TTEST),
         ncol = ncol(InputData))
}
