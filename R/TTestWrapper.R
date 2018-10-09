#' @title T-testing Between Plus and Minus Genes
#'
#' @description .tTest is an internal function used by
#' calculatePhysioMap. It does a t-test between iplus and iminus
#' indexed genes in ReferencesJ.
#'
#' @param ReferencesJ Vector of gene expressions to do statistical test on.
#' @param iplus Index of first group of genes for statistical testing.
#' @param iminus Index of second group of genes for statistical testing.
#' @param STATICResponse Same STATICResponse as in calculatePhysioMap.
#' Check calculatePhysioMap's help for more info.
#'
#' @import stats
#'
#' @return Log2 signed p value of t-test if STATICResponse==FALSE,
#' t-test statisitic if STATICResponse==TRUE.
#'
#' @examples
#'  SimulatedReferenceSpace <-
#'    matrix(
#'      rnorm(n = 100000, mean = 0, sd = 100),
#'      ncol = 10,
#'      dimnames = list(1:10000, 11:20)
#'    )
#'  PhysioSpaceMethods:::.tTest(
#'    ReferencesJ = SimulatedReferenceSpace[, 4],
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

.tTest <- function(ReferencesJ,iplus,iminus,STATICResponse){
    tTestTemp <- t.test(ReferencesJ[iplus], ReferencesJ[iminus])
    if(STATICResponse) {
        Stats <- tTestTemp$statistic
        MaxStat <-
            t.test(sort(ReferencesJ, decreasing = TRUE)[seq_along(iplus)],
                    sort(ReferencesJ)[seq_along(iminus)])$statistic
        return(Stats/MaxStat)
    }
    sgn <- -sign(tTestTemp$statistic)
    pval = tTestTemp$p.value
    return(log2(pval) * sgn)
}
