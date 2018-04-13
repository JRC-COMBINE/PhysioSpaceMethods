#' @title Creates and Uses Physio Spaces as a dimension reduction mapping
#'
#' @description This package uses 'Big Data' to make robust 'Physiological Vectors' in N dimensions spaces, with which you will map new data to extract infromation from a big high dimensional confusing new data.
#'
#' @param ReferencesJ,iplus,iminus,STATICResponse
#'
#' @return NULL
#'
#' @examples tTestWrapper()
#'
#' @export tTestWrapper

#
tTestWrapper <- function(ReferencesJ,iplus,iminus,STATICResponse){
  tTestTemp <- t.test(ReferencesJ[iplus], ReferencesJ[iminus])
  if(STATICResponse) return(tTestTemp$statistic) #ToDo: Want to have values from -1 to 1
  sgn <- -sign(tTestTemp$statistic)
  pval = tTestTemp$p.value
  return(log10(pval) * sgn)
}
