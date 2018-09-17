#' @title Preparing InputData for Calculation
#'
#' @description inptPreparer is an internal function used by
#' inptChecker to set up row- and column names of
#' 'InputData', and convert it to a matrix, if necessary.
#'
#' @inheritParams calculatePhysioMap
#'
#' @importFrom SummarizedExperiment rowData colData assay
#'
#' @return inptPreparer returns InputData as a matrix
#' with proper row- and column names.
#'
#' @examples
#'  SimulatedGeneExpressionData <- matrix(
#'    rnorm(n = 100000, mean = 0,
#'          sd = 100),
#'    ncol = 10,
#'    dimnames = list(1:10000, 1:10)
#'  )
#'  SimulatedGeneExpressionData_checked <-
#'      inptPreparer(SimulatedGeneExpressionData)
#'
#'  library(SummarizedExperiment)
#'  SimulatedGeneExpressionData_SE <- SummarizedExperiment(
#'      assays = list(GEX = SimulatedGeneExpressionData),
#'      rowData = data.frame("EntrezID" =
#'                               rownames(SimulatedGeneExpressionData)),
#'      colData = data.frame("SampleName" =
#'                               colnames(SimulatedGeneExpressionData))
#'  )
#'  SimulatedGeneExpressionData_SE_checked <-
#'      inptPreparer(SimulatedGeneExpressionData_SE)
#'
#' @export
#
inptPreparer <- function(InputData){
    UseMethod("inptPreparer")
}


#' @export
inptPreparer.matrix <- function(InputData){
    if(is.null(rownames(InputData))){
        stop("Rownames of 'InputData' need to be assigned!")
    }
    return(InputData)
}


#' @export
inptPreparer.SummarizedExperiment <- function(InputData){
    InputDataMat <- assay(InputData)
    AvailAnnots <- colData(InputData)
    if(!hasName(AvailAnnots, "SampleName")){
        warning("InputData is missing 'SampleName' ",
                "in its colData. For more info about ",
                "acceptable InputData format check ",
                "the help of 'calculatePhysioMap'")
    }
    colnames(InputDataMat) <- AvailAnnots$SampleName
    AvailGeneAnnots <- rowData(InputData)
    if(!hasName(AvailGeneAnnots, "EntrezID")){
        stop("InputData is missing 'EntrezID' ",
                "in its rowData. For more info about ",
                "acceptable InputData format check ",
                "the help of 'calculatePhysioMap'")
    }
    rownames(InputDataMat) <- AvailGeneAnnots$EntrezID
    return(InputDataMat)
}


#' @export
inptPreparer.default <- function(InputData){
    stop("'InputData' of the class ",
                class(InputData),
                " is not supported!")
}
