#' @title Creates and Uses Physio Spaces as a dimension reduction mapping
#'
#' @description This package uses 'Big Data' to make robust 'Physiological Vectors' in N dimensions spaces, with which you will map new data to extract infromation from a big high dimensional confusing new data.
#'
#' @param InputData,References
#'
#' @return NULL
#'
#' @examples inputChecker()
#'
#' @export inputChecker

#
inputChecker <- function(InputData, References){
  #In case inputs were vectors
  if(!is.matrix(InputData)) stop("'calculatePhysioMap' expects a matrix for InputData")
  if(!is.matrix(References)) stop("'calculatePhysioMap' expects a matrix for References")
  #Imputing missing values:
  if(anyNA(InputData)) InputData <- imputeMissingGeneExpression(InputData)
  if(anyNA(References)) References <- imputeMissingGeneExpression(References)
  #Other checks:
  if (!all(rownames(InputData) == rownames(References))) {
    message("Rows of InputData doesn't match rows of References, trying to match them...")
    matchedIndxes <- match(rownames(InputData), rownames(References))
    commonRowsSum <- sum(!is.na(matchedIndxes))
    if(commonRowsSum < 0.1*min(nrow(InputData), nrow(References))) stop("Less than 10% of rows could be match! aborting...")
    if(commonRowsSum < 200) stop("Less than 200 rows could be match! which is not enough for PhysioSpace, aborting...")
    InputData <- InputData[!is.na(matchedIndxes),]
    References <- References[na.omit(matchedIndxes),]
    message(paste("Matching done, with",nrow(References),"common rows between InputData and References."))
  }
  if (min(InputData, na.rm = T) >= 0)
    warning("You didn't provide relative values in InputData??!")
  #Returning new datasets:
  assign("InputData", value = InputData, envir = parent.frame())
  assign("References", value = References, envir = parent.frame())
}
