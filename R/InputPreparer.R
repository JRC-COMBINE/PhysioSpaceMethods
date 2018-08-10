#' @title Preparing InputData for Calculation
#'
#' @description inputPreparer is an internal function used by
#' inputChecker to set up row- and column names of
#' 'InputData', and convert it to a matrix, if necessary.
#'
#' @inheritParams calculatePhysioMap
#'
#' @importFrom SummarizedExperiment rowData colData assay
#'
#' @return inputPreparer returns InputData as a matrix
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
#'      inputPreparer(SimulatedGeneExpressionData)
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
#'      inputPreparer(SimulatedGeneExpressionData_SE)
#'
#' @export
#
inputPreparer <- function(InputData){
    UseMethod("inputPreparer")
}


#' @export
inputPreparer.matrix <- function(InputData){
    if(is.null(rownames(InputData))){
        stop("Rownames of 'InputData' need to be assigned!")
    }
    return(InputData)
}


#' @export
inputPreparer.SummarizedExperiment <- function(InputData){
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
inputPreparer.default <- function(InputData){
    stop("'InputData' of the class ",
                class(InputData),
                " is not supported!")
}
