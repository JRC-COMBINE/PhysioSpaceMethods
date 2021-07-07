#' @title Imputing Missing Data
#'
#' @description .imputeMissingGeneExpression is an internal function called
#' by inputChecker. It uses different methods to impute any missing value
#' of the InputData or Space (or else missing values would break the
#' pipeline in calculatePhysioMap).
#'
#' @param InptGEX Input matrix with missng values.
#' @param METHOD Method to use in imputation. Available methods are
#' KNN and PCA. Default is 'PCA'.
#'
#' @return A matrix with the same dimensions as InptGEX, with missing
#' values imputed.
#'
#' @import utils
#' @importFrom impute impute.knn
#' @importFrom missMDA imputePCA estim_ncpPCA
#'
#' @examples
#' \dontrun{
#'  MatToImpute <-
#'    matrix(
#'      rnorm(n = 100000, mean = 0, sd = 100),
#'      ncol = 10,
#'      dimnames = list(1:10000, 1:10)
#'    )
#'  MatToImpute[sample(x = 1:length(MatToImpute),
#'                     size = length(MatToImpute) / 20)] <- NA
#'  ImputedMat <-
#'    PhysioSpaceMethods:::.imputeMissingGeneExpression(InptGEX = MatToImpute,
#'    METHOD = "PCA")
#' }
#'
.imputeMissingGeneExpression <- function(InptGEX, METHOD="PCA"){
    if(METHOD=="KNN"){
        ##KNN is shown to be among the worse but popular methods of
        #gene expression imputation, I'll use it for now but have to
        #change later Prob'ly won't work for RNA-seq either -> another
        #reason to change to a better method
        return(impute.knn(data = InptGEX)$data)
    } else if(METHOD=="PCA") {
        res.comp <- imputePCA(InptGEX,ncp = estim_ncpPCA(InptGEX,ncp.max=7)$ncp)
        return(res.comp$completeObs)
    } else {
        stop("Imputation method ", METHOD," is not implemented.")
    }
}
