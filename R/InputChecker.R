#' @title Checking calculatePhysioMap Inputs
#'
#' @description inputChecker is an internal function used by calculatePhysioMap to ckeck the format of inputs, 'InputData' and 'Space' to be
#' exact. It checks to see 1- both 'InputData' and 'Space' are matrices, 2- if 'InputData' or 'Space' contain any NA, and trying to impute
#' using imputeMissingGeneExpression function if they do, and 3- it matches the rows of 'InputData' or 'Space' based on their row names.
#'
#' @param InputData Same InputData as in calculatePhysioMap. Check calculatePhysioMap's help for more info.
#' @param Space Same Space as in calculatePhysioMap. Check calculatePhysioMap's help for more info.
#' @param ImputationMethod Imputation method to use in case of missing values.
#'
#' @return inputChecker returns corrected 'InputData' and 'Space' directly to the environment it was called from
#' (By assigning new matrices to parent.frame()).
#'
#' @examples require(PhysioSpaces)
#' inputChecker(InputData = HS_LUKK_Space[,100:110], Space=HS_LUKK_Space[,1:10])
#'
#' @export inputChecker

#
inputChecker <- function(InputData, Space, ImputationMethod){
  #In case inputs were vectors
  if(!is.matrix(InputData)) stop("'calculatePhysioMap' expects a matrix for InputData")
  if(!is.matrix(Space)) stop("'calculatePhysioMap' expects a matrix for Space")
  #Imputing missing values:
  if(anyNA(InputData)) InputData <- imputeMissingGeneExpression(InputData, METHOD=ImputationMethod)
  if(anyNA(Space)) Space <- imputeMissingGeneExpression(Space, METHOD=ImputationMethod)
  #Other checks:
  if (!identical(rownames(InputData), rownames(Space))) {
    message("Rows of InputData doesn't match rows of Space, trying to match them...")
    matchedIndxes <- match(rownames(InputData), rownames(Space))
    commonRowsSum <- sum(!is.na(matchedIndxes))
    if(commonRowsSum < 0.1*min(nrow(InputData), nrow(Space))) stop("Less than 10% of rows could be match! aborting...")
    if(commonRowsSum < 200) stop("Less than 200 rows could be match! which is not enough for PhysioSpace, aborting...")
    InputData <- InputData[!is.na(matchedIndxes),]
    Space <- Space[na.omit(matchedIndxes),]
    message(paste("Matching done, with",nrow(Space),"common rows between InputData and Space."))
  }
  if (min(InputData, na.rm = T) >= 0)
    warning("You didn't provide relative values in InputData??!")
  #Returning new datasets:
  assign("InputData", value = InputData, envir = parent.frame())
  assign("Space", value = Space, envir = parent.frame())
}
