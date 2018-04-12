#' @title Creates PhysioSpaces
#'
#' @description This package uses 'Big Data' to make robust 'Physiological Vectors' in N dimensions spaces, with which you will map new data to extract infromation from a big high dimensional confusing new data.
#'
#' @param InputData
#'
#' @return NULL
#'
#' @examples  spaceMaker()
#'
#' @export spaceMaker
#'
spaceMaker <- function(inputMatrix, annotationFrame, DESIGN, CONTRAST, RETURNMODEL = F, MODEL = NULL, RETURNFD = F){
  if(is.null(MODEL)) {
    if(!is.matrix(inputMatrix)) stop("'inputMatrix' should be a matrix!")
    if(!is.data.frame(annotationFrame)) stop("'annotationFrame' should be a data.frame!")
    require(DESeq2)
    MODELStructured <- DESeqDataSetFromMatrix(countData = inputMatrix, colData = annotationFrame,
                                              design = DESIGN)
    MODELCalculated <- DESeq(MODELStructured)
  } else {
    MODELCalculated <- MODEL
  }
  MODELResults <- results(object = MODELCalculated, contrast = CONTRAST)
  if(RETURNMODEL) return(MODELCalculated)
  if(RETURNFD) return(MODELResults$log2FoldChange)
  return((1/MODELResults$pvalue)*sign(MODELResults$log2FoldChange))
}
