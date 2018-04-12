#' @title Creates and Uses Physio Spaces as a dimension reduction mapping
#'
#' @description This package uses 'Big Data' to make robust 'Physiological Vectors' in N dimensions spaces, with which you will map new data to extract infromation from a big high dimensional confusing new data.
#'
#' @param InptGEX
#'
#' @return NULL
#'
#' @examples  imputeMissingGeneExpression()
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
