#' @title Mapping new data into a physiological-space
#'
#' @description calculatePhysioMap computes mapped values of each input sample inside of a space,
#' calculated prior from a compendium of known samples.
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
#' It is not recomemnded to use PARALLEL=TRUE on small datasets since due to large overhead, it could take more time than
#' using PARALLEL=FALSE.
#' @param NumbrOfCores Number of cores to be used when 'PARALLEL' is TRUE. Default is NA which will result in the program using
#' 'all available cores - 1'. Assigning a number higher that the value parallel::detectCores() returns will result in an error.
#' @param TTEST Logical value indicating if t.test should be done in place of the default wilcoxon rank-sum test (more info can be found
#' in the original PhysioSpace: Lenz et. al., PLOS One 2013). Using t.test will speed up calculations. Default value is FALSE.
#' @param STATICResponse Logical value indicating if 'statistic' should be returned rather than the default 'signed p value'.
#' Default value is FALSE.
#' @param ImputationMethod Imputation method to use in case of missing values.
#' @param ParallelMethod Parallel method to use when PARALLEL=TRUE. Two methods are implemented so far: "parCapply"
#' which uses parCapply function of parallel package, and "foreach" which uses foreach package to process in parallel.
#' Speed-wise, foreach has an edge on parCapply, but the latter is more stable. Hence we recommend parCapply.
#' The default value for ParallelMethod is "parCapply".
#'
#' @import progress
#'
#' @details ToDo
#'
#' @return Matrix of mapped 'InputData' values in 'Space', with rows corrisponding to axises of 'Space' and columns representing
#' samples in 'InputData'. Mapped values are signed p value when STATICResponse==FALSE, and are 'statistic' value when
#' STATICResponse==TRUE (more info can be found in the original PhysioSpace: Lenz et. al., PLOS One 2013).
#'
#'
#' @examples
#' SimulatedGeneExpressionData <- matrix(rnorm(n = 100000, mean = 0, sd = 100),ncol = 10, dimnames = list(1:10000,1:10))
#' SimulatedReferenceSpace <- matrix(rnorm(n = 100000, mean = 0, sd = 100),ncol = 10, dimnames = list(1:10000,11:20))
#' calculatePhysioMap(InputData = SimulatedGeneExpressionData, SimulatedReferenceSpace)
#' if(parallel::detectCores()>1){ #More than one core is needed for parallel processing
#' calculatePhysioMap(InputData = SimulatedGeneExpressionData, SimulatedReferenceSpace,
#' PARALLEL=TRUE, NumbrOfCores=2, GenesRatio = 0.01, STATICResponse = FALSE, TTEST = TRUE)
#' }
#'
#' @export
#Pre-function for dispaching the calculating to the write function:
calculatePhysioMap <- function(InputData, Space, GenesRatio = 0.05,
                               PARALLEL = FALSE, NumbrOfCores=NA,
                               TTEST = FALSE, STATICResponse = FALSE,
                               ImputationMethod = "PCA", ParallelMethod = "parCapply"){
  UseMethod("calculatePhysioMap", object = Space)
}


##The main function for calculating PhysioScores:
#' @export
calculatePhysioMap.default <- function(InputData, Space, GenesRatio = 0.05,
                                       PARALLEL = FALSE, NumbrOfCores=NA,
                                       TTEST = FALSE, STATICResponse = FALSE,
                                       ImputationMethod = "PCA", ParallelMethod = "parCapply"){
  #Initializing:
  inputChecker(InputData, Space, ImputationMethod)
  NGenes <- nrow(InputData)
  NSamples <- ncol(InputData)
  physioMap <- matrix(NA, ncol(Space), NSamples)
  if(PARALLEL) cl <- parallelInitializer(NumbrOfCores=NumbrOfCores)
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent eta: :eta",
                         total = ifelse(PARALLEL,NSamples/length(cl),NSamples), clear = FALSE)
  #
  #Main:
  if(PARALLEL) {
    if(ParallelMethod == "parCapply"){
      physioMap <- calculatePhysioMapCore_parCapply(InputData, Space, NSamples,
                                                    GenesRatio, NGenes, STATICResponse, pb, TTEST, cl)
    } else if(ParallelMethod == "foreach"){
      physioMap <- calculatePhysioMapCore_foreach(InputData, Space,
                                                  NSamples, GenesRatio, NGenes, STATICResponse, pb, TTEST)
    } else {
      stop("The specified ParallelMethod is not implemented!")
    }
    stopCluster(cl)
  } else {
    physioMap <- apply(X = InputData,
                       MARGIN = 2, FUN = singleThreadOfPhysioCalc, Space=Space,
                       NSamples=NSamples, GenesRatio=GenesRatio, NGenes=NGenes,
                       STATICResponse=STATICResponse, pb=pb, TTEST=TTEST)
  }

  rownames(physioMap) = colnames(Space)
  colnames(physioMap) = colnames(InputData)
  return(physioMap)
}


##PhysioSpace for gene lists as a Space:
#' @export
calculatePhysioMap.list <- function(InputData, Space, GenesRatio = 0.05,
                                    PARALLEL = FALSE, NumbrOfCores=NA,
                                    TTEST = FALSE, STATICResponse = FALSE,
                                    ImputationMethod = "PCA", ParallelMethod = "parCapply"){

  #Check the format of the input, and correct:(How can I move this function to its own .R file??)
  inputChecker <- function(InputData, Space){
    #In case inputs were vectors
    if(!is.matrix(InputData)) InputData <<- as.matrix(InputData)
    #Imputing missing values:
    if(anyNA(InputData)) InputData <- imputeMissingGeneExpression(InputData)
    #Other checks:
    if(is.null(names(Space))) stop("list-Space should have 'names'!")
    if (min(InputData, na.rm = TRUE) >= 0)
      warning("You didn't provide relative values in InputData??!")
  }

  #Initializing:
  warning("Not implemented yet!")
  return(NULL)
}
