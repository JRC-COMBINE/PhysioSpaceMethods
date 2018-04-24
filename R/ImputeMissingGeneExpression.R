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
#' @import missMDA mice
#'
#' @examples require(PhysioSpaces)
#' MatToImpute <- HS_LUKK_Space[,100:110]
#' MatToImpute[sample(x = 1:length(MatToImpute), size = length(MatToImpute)/20)] <- NA
#' ImputedMat <- imputeMissingGeneExpression(InptGEX = MatToImpute, METHOD = "PCA")
#' ImputedMat2 <- imputeMissingGeneExpression(InptGEX = MatToImpute, METHOD = "KNN")
#' ImputedMat3 <- imputeMissingGeneExpression(InptGEX = MatToImpute, METHOD = "Regression")
#'
#' @export imputeMissingGeneExpression

imputeMissingGeneExpression <- function(InptGEX, METHOD="KNN"){
  if(METHOD=="KNN"){
    ##KNN is shown to be among the worse but popular methods of gene expression imputation, I'll use it for now but have to change later
    #Prob'ly won't work for RNA-seq either -> another reason to change to a better method
    loadedSuccessFully <- require(impute, quietly = T)
    if(!loadedSuccessFully){
      warning("Impute package wasn't found, will try to install it from Bioconductor...")
      source("https://bioconductor.org/biocLite.R")
      biocLite("impute")
      library(impute)
    }
    invisible(capture.output(imputedInpt <- impute::impute.knn(data = InptGEX)$data)) #'invisible+capture.output' so I could suppress all irritating cat messages of impute.knn
    return(imputedInpt)
  } else if(METHOD=="PCA") {
    res.comp <- imputePCA(InptGEX,ncp = estim_ncpPCA(InptGEX,ncp.max=7)$ncp)
    return(res.comp$completeObs)
  } else if(METHOD=="Regression") {
    suppressPackageStartupMessages(require(mice)) #Although mice is loaded via namespace cause it's part of
    # PhysioSpaceMethods imports, I have to manually load it or else next line may break (weird implementation in mice)
    # Also have to convert to data.frame and take care of dimnames so mice wouldn't break:
    dimNameKeeps <- dimnames(InptGEX)
    InptGEX <- as.data.frame(InptGEX, row.names = make.names(rownames(InptGEX), unique = T),
                             col.names = make.names(colnames(InptGEX), unique = T))
    imp <- mice(InptGEX, method = "norm.predict", m = 1, printFlag=F)
    imputed <- as.matrix(complete(imp))
    dimnames(imputed) <- dimNameKeeps
    return(imputed)
  } else {
    stop(paste("Imputation method",METHOD,"not implemented!"))
  }
}
