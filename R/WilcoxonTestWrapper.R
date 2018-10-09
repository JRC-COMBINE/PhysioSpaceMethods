#' @title Wilcoxon Rank Sum testing Between Plus and Minus Genes
#'
#' @description .wilTest is an internal function used by
#' calculatePhysioMap that does rank sum test
#' (equivalent to Mann-Whitney test) between iplus and
#' iminus indexed genes in ReferencesJ.
#'
#' @param ReferencesJ Vector of gene expressions to do statistical test on.
#' @param iplus Index of first group of genes for statistical testing.
#' @param iminus Index of second group of genes for statistical testing.
#' @param STATICResponse Same STATICResponse as in calculatePhysioMap.
#' Check calculatePhysioMap's help for more info.
#'
#' @import stats
#'
#' @return Log2 signed p value of Rank sum test if STATICResponse==FALSE,
#' Rank sum statisitic normalized between -1 and 1
#' if STATICResponse==TRUE.
#'
#' @examples
#'  SimulatedReferenceSpace <-
#'    matrix(
#'      rnorm(n = 100000, mean = 0, sd = 100),
#'      ncol = 10,
#'      dimnames = list(1:10000, 11:20)
#'    )
#'  PhysioSpaceMethods:::.wilTest(
#'    ReferencesJ = SimulatedReferenceSpace[, 9],
#'    iplus = sample(
#'      1:nrow(SimulatedReferenceSpace),
#'      size = nrow(SimulatedReferenceSpace) / 20
#'    ),
#'    iminus = sample(
#'      1:nrow(SimulatedReferenceSpace),
#'      size = nrow(SimulatedReferenceSpace) / 20
#'    ),
#'    STATICResponse = FALSE
#'  )
#'

.wilTest <- function(ReferencesJ,iplus,iminus,STATICResponse){
    wilTestTemp <- wilcox.test(ReferencesJ[iplus], ReferencesJ[iminus])
    if(STATICResponse) {
        return(2*((wilTestTemp$statistic/(length(iplus)*length(iminus)))-0.5))
        #Want to have values from -1 to 1
    }
    sgn <- sign(-(wilTestTemp$statistic - (length(iplus)*length(iminus)) / 2))
    pval = wilTestTemp$p.value
    return(log2(pval) * sgn)
}
