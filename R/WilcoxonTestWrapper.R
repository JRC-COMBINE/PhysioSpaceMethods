#' @title Creates and Uses Physio Spaces as a dimension reduction mapping
#'
#' @description This package uses 'Big Data' to make robust 'Physiological Vectors' in N dimensions spaces, with which you will map new data to extract infromation from a big high dimensional confusing new data.
#'
#' @param ReferencesJ,iplus,iminus,WilxSTATICResponse
#'
#' @return NULL
#'
#' @examples wilTestWrapper()
#'
#' @export wilTestWrapper

#
wilTestWrapper <- function(ReferencesJ,iplus,iminus,WilxSTATICResponse){
  wilTestTemp <- wilcox.test(ReferencesJ[iplus], ReferencesJ[iminus])
  if(WilxSTATICResponse) return(2*((wilTestTemp$statistic/(length(iplus)*length(iminus)))-0.5)) #Want to have values from -1 to 1
  sgn <- sign(-(wilTestTemp$statistic - (length(iplus)*length(iminus)) / 2))
  pval = wilTestTemp$p.value
  return(log10(pval) * sgn)
}
