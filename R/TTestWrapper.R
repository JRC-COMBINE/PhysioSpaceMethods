#' @title Creates and Uses Physio Spaces as a dimension reduction mapping
#'
#' @description This package uses 'Big Data' to make robust 'Physiological Vectors' in N dimensions spaces, with which you will map new data to extract infromation from a big high dimensional confusing new data.
#'
#' @param ReferencesJ,iplus,iminus
#'
#' @return NULL
#'
#' @examples TTestWrapper()
#'
#' @export TTestWrapper

#
TTestWrapper <- function(ReferencesJ,iplus,iminus,...){ #... is needed so apply with WilxSTATICResponse input won't break this function
  TTestTemp <- t.test(ReferencesJ[iplus], ReferencesJ[iminus])
  sgn <- -sign(TTestTemp$statistic)
  pval = TTestTemp$p.value
  return(log10(pval) * sgn)
}
