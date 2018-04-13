#' @title Mapping new data into a physiological-space
#'
#' @description calculatePhysioMap computes mapped values of each input sample inside of a space,
#' calculated prior from a compendium of known samples (Check ?SpaceMaker for more info).
#'
#' @param InputData A matrix of input gene expressions to be analyzed, with genes as rows and samples as columns. Corresponding Entrez Gene
#' IDs must be assigned to 'rownames' of the matrix, and name of each sample/column should be written in 'colnames'. REMEMBER that the gene
#' expressions in 'InputData' should be relative; e.g. fold change or signed p value of a statistical test.
#' @param Space The space in which the 'InputData' will be mapped. Just as 'InputData',
#' it should be a matrix with genes as rows and samples as columns, with corresponding Entrez Gene
#' IDs in 'rownames' of the matrix, and name of each axis of the space written in 'colnames'.
#' @param GenesRatio The ratio of gene expression values to be considered in the calculation. In the calculations only the highest and lowest
#' GenesRatio*100 percent of values in each sample are used, since signal to noise ratio in gene expression values
#' has a direct relation to the relative magnitute of expressions and out noisy inputs need to be filtered out.
#' GenesRatio should be a numerical value between 0 and 1. Default value is 0.05.
#' @param PARALLEL Logical value indicating if calculations should be done in parallel. Default value is FALSE.
#' @param NumbrOfCores Number of cores to be used when 'PARALLEL' is TRUE. Default is NA which will result in the program using
#' 'all available cores - 1'.
#' @param TTEST Logical value indicating if t.test should be done in place of the default wilcoxon rank-sum test (more info can be found
#' in the original PhysioSpace: Lenz et. al., PLOS One 2013). Using t.test will speed up calculations. Default value is FALSE.
#' @param STATICResponse Logical value indicating if 'statistic' should be returned rather than the default 'signed p value'.
#' Default value is FALSE.
#'
#' @details ToDo
#'
#' @return Matrix of mapped 'InputData' values in 'Space', with rows corrisponding to axises of 'Space' and columns representing
#' samples in 'InputData'. Mapped values are signed p value when STATICResponse==FALSE, and are 'statistic' value when
#' STATICResponse==TRUE (more info can be found in the original PhysioSpace: Lenz et. al., PLOS One 2013).
#'
#'
#' @examples require(PhysioSpaces)
#' calculatePhysioMap(InputData = HS_LUKK_Space[,100:110], Space=HS_LUKK_Space[,1:10])
#' calculatePhysioMap(InputData = HS_LUKK_Space[,100:110], Space=HS_U133Plus2_Space[,1:10],
#' PARALLEL=T, NumbrOfCores=4, GenesRatio = 0.01, STATICResponse = F, TTEST = T)
#'
#' @export calculatePhysioMap
#Pre-function for dispaching the calculating to the write function:
calculatePhysioMap <- function(InputData, Space, GenesRatio = 0.05,
                               PARALLEL = FALSE, NumbrOfCores=NA,
                               TTEST = FALSE, STATICResponse = F){
  UseMethod("calculatePhysioMap", object = Space)
}

#' @export calculatePhysioMap.default
##The main function for calculating PhysioScores:
calculatePhysioMap.default <- function(InputData, Space, GenesRatio = 0.05,
                                       PARALLEL = FALSE, NumbrOfCores=NA,
                                       TTEST = FALSE, STATICResponse = F){
  #Initializing:
  inputChecker(InputData, Space)
  NGenes <- nrow(InputData)
  NSamples <- ncol(InputData)
  physioMap <- matrix(NA, ncol(Space), NSamples)
  suppressPackageStartupMessages(require(progress))
  suppressPackageStartupMessages(require(doParallel))
  if(PARALLEL) cl <- parallelInitializer(NumbrOfCores=NumbrOfCores)
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent eta: :eta",
                         total = ifelse(PARALLEL,NSamples/length(cl),NSamples), clear = FALSE)
  #
  #Main:
  suppressWarnings( #So dopar won't make a warning if PARALELL=F and cl doesn't exist
  physioMap <-
    foreach(i = 1:NSamples, .combine='cbind', .final=as.matrix, .export=c("tTestWrapper","wilTestWrapper"), .packages = "PhysioSpaceMethods") %dopar% {
      tempDiff <- InputData[, i]
      if (!is.null(GenesRatio)) {
        numGenes = round(NGenes * GenesRatio)
      } else{
        stop("Don't have Number of genes to compare??!?!")
      }
      ordDiff = order(tempDiff)
      if (!is.null(numGenes)) {
        iplus = ordDiff[(NGenes - numGenes + 1):NGenes]
        iminus = ordDiff[1:numGenes]
      }
      pb$tick()
      apply(X = Space, MARGIN = 2, FUN = if(TTEST) tTestWrapper else wilTestWrapper,
            iplus=iplus, iminus = iminus, STATICResponse = STATICResponse) # this apply can be written as simple for loop also, results
      ## seems to be the same, but writing a nested foreach loop may speed up stuff -> will have to try later
    }
  )
  if(PARALLEL) stopCluster(cl)

  rownames(physioMap) = colnames(Space)
  colnames(physioMap) = colnames(InputData)
  return(physioMap)
}

#' @export calculatePhysioMap.list
##PhysioSpace for gene lists as a Space:
calculatePhysioMap.list <- function(InputData, Space, GenesRatio = 0.05,
                                    PARALLEL = FALSE, NumbrOfCores=NA,
                                    TTEST = FALSE, STATICResponse = F){

  #Check the format of the input, and correct:(How can I move this function to its own .R file??)
  inputChecker <- function(InputData, Space){
    #In case inputs were vectors
    if(!is.matrix(InputData)) InputData <<- as.matrix(InputData)
    #Imputing missing values:
    if(anyNA(InputData)) InputData <- imputeMissingGeneExpression(InputData)
    #Other checks:
    if(is.null(names(Space))) stop("list-Space should have 'names'!")
    if (min(InputData, na.rm = T) >= 0)
      warning("You didn't provide relative values in InputData??!")
  }

  #Initializing:
  warning("Not implemented yet!")
  return(NULL)
}
