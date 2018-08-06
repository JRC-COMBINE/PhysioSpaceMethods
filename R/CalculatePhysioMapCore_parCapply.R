#' @title Mapping new data into a physiological-space using parCapply of parallel package
#'
#' @description calculatePhysioMapCore_parCapply is an internal function of calculatePhysioMap,
#' computing the main mapping using pacCapply of parallel package. We don't recommend the use of
#' calculatePhysioMapCore_parCapply outside of calculatePhysioMap.
#'
#' @import progress parallel
#'
#' @param InputData A matrix of input gene expressions to be analyzed, with genes as rows and samples as columns.
#' Corresponding Entrez Gene
#' IDs must be assigned to 'rownames' of the matrix, and name of each sample/column should be written in 'colnames'.
#' REMEMBER that the gene
#' expressions in 'InputData' should be relative; e.g. fold change or signed p value of a statistical test.
#' @param Space The space in which the 'InputData' will be mapped. Just as 'InputData',
#' it should be a matrix with genes as rows and samples as columns, with corresponding Entrez Gene
#' IDs in 'rownames' of the matrix, and name of each axis of the space written in 'colnames'.
#' @param GenesRatio The ratio of gene expression values to be considered in the calculation. In high dimensional omics data,
#' signal to noise ratio has a direct relation with the relative magnitude of expressions. We aim to remove the noisy genes,
#' hence we only keep the "GenesRatio*100" percent highest and lowest gene expression values of each sample.
#' GenesRatio should be a numerical value between 0 and 1. Default value is 0.05.
#' @param TTEST Logical value indicating if t.test should be done in place of the default wilcoxon rank-sum test (more info can be found
#' in the original PhysioSpace: Lenz et. al., PLOS One 2013). Using t.test will speed up calculations. Default value is FALSE.
#' @param STATICResponse Logical value indicating if 'statistic' should be returned rather than the default 'signed p value'.
#' Default value is FALSE.
#' @param NSamples Number of samples (columns) in InputData.
#' @param NGenes Number of genes (rows) in InputData or Space.
#' @param pb Progress bar made by progress::progress_bar$new.
#' @param cl Cluster made by parallel::makeCluster.
#'
#' @return  Matrix of mapped 'InputData' values in 'Space', with rows corresponding to axes of 'Space' and columns representing
#' samples in 'InputData'. Mapped values are signed p value when STATICResponse==FALSE, and are 'statistic' value when
#' STATICResponse==TRUE (more info can be found in the original PhysioSpace paper: Lenz et. al., PLOS One 2013).
#'
# #' @export #Not exporting since it's an internal function
calculatePhysioMapCore_parCapply <- function(InputData, Space, NSamples, GenesRatio, NGenes, STATICResponse, pb, TTEST, cl){
  matrix(parCapply(cl = cl, x = InputData, FUN = singleThreadOfPhysioCalc, Space=Space,
                   GenesRatio=GenesRatio, NGenes=NGenes,
                   STATICResponse=STATICResponse, pb=pb, TTEST=TTEST),
         ncol = ncol(InputData))
}
