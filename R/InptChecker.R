#' @title Checking calculatePhysioMap Inputs
#'
#' @description inptChecker is an internal function used by
#' calculatePhysioMap to check the format of inputs, 'InputData'
#' and 'Space' to be exact. It checks to see 1- both 'InputData'
#' and 'Space' are matrices, 2- if 'InputData' or 'Space' contain
#' any NA, and trying to impute using imputeMissingGeneExpression
#' function if they do, and 3- it matches the rows of 'InputData'
#' or 'Space' based on their row names.
#'
#' @inheritParams calculatePhysioMap
#'
#' @import stats
#'
#' @return inptChecker returns corrected 'InputData' and 'Space'
#' directly to the environment it was called from
#' (By assigning new matrices to parent.frame()).
#'
#' @examples
#'  SimulatedGeneExpressionData <-
#'    matrix(
#'      rnorm(n = 100000, mean = 0, sd = 100),
#'      ncol = 10,
#'      dimnames = list(1:10000, 1:10)
#'    )
#'  inptChecker(InputData = SimulatedGeneExpressionData[, 1:5],
#'               Space = SimulatedGeneExpressionData[sample(1:10000), 6:10])
#'
#' @export
#
inptChecker <- function(InputData, Space, ImputationMethod){
    #Space Checking:
    if(!is.matrix(Space))
        stop("'calculatePhysioMap' expects a matrix for Space"
        )
    #Imputing missing values:
    if(anyNA(InputData)) InputData <-
            imputeMissingGeneExpression(InputData,
                                        METHOD=ImputationMethod)
    if(anyNA(Space)) Space <-
            imputeMissingGeneExpression(Space,
                                        METHOD=ImputationMethod)
    #Setting up and preparing InputData as a Matrix (if needed):
    InputData <- inptPreparer(InputData)
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
        InputData <- InputData[!is.na(matchedIndxes),]
        Space <- Space[na.omit(matchedIndxes),]
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
