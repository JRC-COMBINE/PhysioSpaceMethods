#' @title Imputing Missing Data
#'
#' @description imputeMissingGeneExpression is an internal function called by inputChecker. It uses KNN method to impute
#' any missing value of the InputData or Space (or else missing values would break the pipeline in calculatePhysioMap).
#'
#' KNN is shown to be not the best method to impute gene expression data, so we recommend that user imputes any missing
#' value themselves before using calculatePhysioMap.
#'
#' @param InptGEX Input matrix with missng values
#' @param METHOD Method to use in imputation. For now only KNN is available.
#'
#' @return A matrix with the same dimensions as InptGEX, with missing values imputed.
#'
#' @examples require(PhysioSpaces)
#' MatToImpute <- HS_LUKK_Space[,100:110]
#' MatToImpute[sample(x = 1:length(MatToImpute), size = length(MatToImpute)/20)] <- NA
#' ImputedMat <- imputeMissingGeneExpression(InptGEX = MatToImpute, METHOD = "KNN")
#'
#' @export imputeMissingGeneExpression

imputeMissingGeneExpression <- function(InptGEX, METHOD="KNN"){
  if(METHOD=="KNN"){
    ##KNN is shown to be among the worse but popular methods of gene expression imputation, I'll use it for now but have to change later
    #Prob'ly won't work for RNA-seq either -> another reason to change to a better method
    invisible(capture.output(imputedInpt <- impute::impute.knn(data = InptGEX)$data)) #'invisible+capture.output' so I could suppress all irritating cat messages of impute.knn
    return(imputedInpt)
  } else {
    stop(paste("Imputation method",METHOD,"not implemented yet"))
  }
}
