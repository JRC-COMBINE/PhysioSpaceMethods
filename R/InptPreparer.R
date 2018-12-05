#' @title Preparing InputData for Calculation
#'
#' @description .inptPreparer is an internal function used by
#' .inptChecker and spaceMaker to set up row and column names of
#' 'InputData', and convert it to a matrix, if necessary.
#'
#' @param InputData A matrix, SummarizedExperiment object or a list,
#' based on the gene expression data (or any other type of high
#' dimensional data, e.g. protein abundance, SNP, Methylation, etc.),
#' to be analysed. InputData has to have a specific format to be
#' properly analysed, these requirements are thoroughly explained
#' in the 'Details' section of calculatePhysioMap() function.
#'
#' @importFrom SummarizedExperiment rowData colData assay
#'
#' @return .inptPreparer returns InputData as a matrix
#' with proper row and column names.
#'
#' @examples
#' \dontrun{
#'  SimulatedGeneExpressionData <- matrix(
#'    rnorm(n = 100000, mean = 0,
#'          sd = 100),
#'    ncol = 10,
#'    dimnames = list(1:10000, 1:10)
#'  )
#'  SimulatedGeneExpressionData_checked <-
#'      PhysioSpaceMethods:::.inptPreparer(SimulatedGeneExpressionData)
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
#'      PhysioSpaceMethods:::.inptPreparer(SimulatedGeneExpressionData_SE)
#' }
#'
.inptPreparer <- function(InputData){
    UseMethod(".inptPreparer")
}

#'@export
.inptPreparer.matrix <- function(InputData){
    if(is.null(rownames(InputData))){
        stop("Rownames of 'InputData' need to be assigned!")
    }
    return(InputData)
}

#'@export
.inptPreparer.SummarizedExperiment <- function(InputData){
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

#'@export
.inptPreparer.default <- function(InputData){
    stop("'InputData' of the class ",
                class(InputData),
                " is not supported!")
}
