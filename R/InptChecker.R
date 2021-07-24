#' @title Checking calculatePhysioMap Inputs
#'
#' @description .inptChecker is an internal function used by
#' calculatePhysioMap for checking the format of inputs, 'InputData'
#' and 'Space' to be exact. It 1- checks to see if both 'InputData'
#' and 'Space' are matrices, and 2- matches the rows of 'InputData'
#' or 'Space' based on their row names.
#'
#' @param InputData A matrix, SummarizedExperiment object or a list,
#' based on the gene expression data (or any other type of high
#' dimensional data, e.g. protein abundance, SNP, Methylation, etc.),
#' to be analysed. InputData has to have a specific format to be
#' properly analysed, these requirements are thoroughly explained
#' in the 'Details' section of calculatePhysioMap() function.
#'
#' @inheritParams calculatePhysioMap
#'
#' @import stats
#'
#' @return .inptChecker returns corrected 'InputData' and 'Space'
#' directly to the environment it was called from
#' (By assigning new matrices to parent.frame()).
#'
#' @examples
#' \dontrun{
#'  SimulatedGeneExpressionData <-
#'    matrix(
#'      rnorm(n = 100000, mean = 0, sd = 100),
#'      ncol = 10,
#'      dimnames = list(1:10000, 1:10)
#'    )
#'  PhysioSpaceMethods:::.inptChecker(InputData =
#'                                          SimulatedGeneExpressionData[, 1:5],
#'               Space = SimulatedGeneExpressionData[sample(1:10000), 6:10])
#' }
#'
.inptChecker <- function(InputData, Space){
    #Space Checking:
    if(!is.matrix(Space))
        stop("'calculatePhysioMap' expects a matrix for Space"
        )
    #Setting up and preparing InputData as a Matrix (if needed):
    InputData <- .inptPreparer(InputData)
    #Other checks:
    if (!identical(rownames(InputData), rownames(Space))) {
        message("Rows of InputData doesn't match rows of Space,",
                " trying to match them...")
        matchedIndxes <- match(rownames(InputData), rownames(Space))
        commonRowsSum <- sum(!is.na(matchedIndxes))
        if(commonRowsSum < 0.1*min(nrow(InputData), nrow(Space))) {
            stop("Less than 10% of rows could be matched! aborting...")
        }
        if(commonRowsSum < 200) {
            stop(
                "Less than 200 rows could be matched!",
                " which is not enough for PhysioSpace, aborting..."
            )
        }
        InputData <- InputData[!is.na(matchedIndxes),,drop=FALSE]
        Space <- Space[na.omit(matchedIndxes),,drop=FALSE]
        message("Matching done, with ", nrow(Space),
                        " common rows between InputData and Space.")
    }
    if (min(InputData, na.rm = TRUE) >= 0)
        warning("'calculatePhysioMap' expects relative values in",
                " InputData, please make sure values in 'InputData'",
                " are actually relative values.")
    #Returning new datasets:
    assign("InputData", value = InputData, envir = parent.frame())
    assign("Space", value = Space, envir = parent.frame())
}
