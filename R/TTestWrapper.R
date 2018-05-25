#' @title T-testing Between Plus and Minus Genes
#'
#' @description tTestWrapper is an internal function used by calculatePhysioMap that does t-test between iplus and iminus indexed genes in
#' ReferencesJ.
#'
#' @param ReferencesJ Vector of gene expressions to do statistical test on.
#' @param iplus Index of first group of genes for statistical testing.
#' @param iminus Index of first group of genes for statistical testing.
#' @param STATICResponse Same STATICResponse as in calculatePhysioMap. Check calculatePhysioMap's help for more info.
#'
#' @return Log2 signed p value of t-test if STATICResponse==FALSE, t-test statisitic if STATICResponse==TRUE.
#'
#' @examples require(PhysioSpaces)
#' tTestWrapper(ReferencesJ = HS_LUKK_Space[,91], iplus = sample(1:nrow(HS_LUKK_Space), size = nrow(HS_LUKK_Space)/20),
#' iminus = sample(1:nrow(HS_LUKK_Space), size = nrow(HS_LUKK_Space)/20), STATICResponse = F)
#'
#' @export tTestWrapper

#
tTestWrapper <- function(ReferencesJ,iplus,iminus,STATICResponse){
  tTestTemp <- t.test(ReferencesJ[iplus], ReferencesJ[iminus])
  if(STATICResponse) {
    Stats <- tTestTemp$statistic
    MaxStat <- t.test(sort(ReferencesJ,decreasing = T)[1:length(iplus)], sort(ReferencesJ)[1:length(iminus)])$statistic
    return(Stats/MaxStat)
  }
  sgn <- -sign(tTestTemp$statistic)
  pval = tTestTemp$p.value
  return(log10(pval) * sgn)
}
