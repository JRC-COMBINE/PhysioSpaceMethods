#' @title Dispatching the parallel calculation among core functions
#'
#' @description calculatePhysioMapCores is an internal
#' function of calculatePhysioMap, that calls the proper internal
#' function for computing the main mapping in parallel. We don't
#' recommend the use of calculatePhysioMapCores outside of
#' calculatePhysioMap.
#'
#'
#' @inheritParams calculatePhysioMapCorePC
#'
#' @param ParallelMethod Parallel method to use.
#' Two methods are implemented so far: "parCapply"
#' which uses parCapply function of parallel package, and "foreach" which
#' uses foreach package to process in parallel.
#' Speed-wise, foreach has an edge on parCapply, but the latter
#' is more stable. Hence, we recommend to use parCapply.
#' The default value for ParallelMethod is "parCapply".
#'
#' @return  Matrix of mapped 'InputData' values in 'Space', with
#' rows corresponding to axes of 'Space' and columns representing
#' samples in 'InputData'.
#'
# #' @export #Not exporting since it's an internal function
calculatePhysioMapCores <- function(InputData, Space,
                                                NSamples, GenesRatio,
                                                NGenes, STATICResponse,
                                                pb, TTEST, cl,
                                                ParallelMethod="parCapply"){
    switch(ParallelMethod,
            parCapply = calculatePhysioMapCorePC(InputData = InputData,
                                                Space=Space,
                                                NSamples=NSamples,
                                                GenesRatio=GenesRatio,
                                                NGenes=NGenes,
                                                STATICResponse=STATICResponse,
                                                pb=pb, TTEST=TTEST, cl=cl),
            foreach = calculatePhysioMapCoreFE(InputData = InputData,
                                                Space=Space,
                                                NSamples=NSamples,
                                                GenesRatio=GenesRatio,
                                                NGenes=NGenes,
                                                STATICResponse=STATICResponse,
                                                pb=pb, TTEST=TTEST),
            stop("The specified ParallelMethod is not implemented."))
}
