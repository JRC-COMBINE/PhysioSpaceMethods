#' @title Imputing Missing Data
#'
#' @description imputeMissingGeneExpression is an internal function called by inputChecker. It uses different methods to impute
#' any missing value of the InputData or Space (or else missing values would break the pipeline in calculatePhysioMap).
#'
#' @param InptGEX Input matrix with missng values
#' @param METHOD Method to use in imputation.
#'
#' @return A matrix with the same dimensions as InptGEX, with missing values imputed.
#'
#' @import missMDA DMwR utils
#'
#' @examples
#'  MatToImpute <-
#'    matrix(
#'      rnorm(n = 100000, mean = 0, sd = 100),
#'      ncol = 10,
#'      dimnames = list(1:10000, 1:10)
#'    )
#'  MatToImpute[sample(x = 1:length(MatToImpute),
#'                     size = length(MatToImpute) / 20)] <- NA
#'  ImputedMat <-
#'    imputeMissingGeneExpression(InptGEX = MatToImpute, METHOD = "PCA")
#'  ImputedMat2 <-
#'    imputeMissingGeneExpression(InptGEX = MatToImpute, METHOD = "KNN")
#'
#' @export imputeMissingGeneExpression

imputeMissingGeneExpression <- function(InptGEX, METHOD="PCA"){
  if(METHOD=="KNN"){
    ##KNN is shown to be among the worse but popular methods of gene expression imputation, I'll use it for now but have to change later
    #Prob'ly won't work for RNA-seq either -> another reason to change to a better method
    return(knnImputation(data = InptGEX))
  } else if(METHOD=="PCA") {
    res.comp <- imputePCA(InptGEX,ncp = estim_ncpPCA(InptGEX,ncp.max=7)$ncp)
    return(res.comp$completeObs)
  } else {
    stop(paste("Imputation method",METHOD,"not implemented!"))
  }
}
