#' @title Mapping new data into a physiological-space
#'
#' @description calculatePhysioMap computes mapped values of each input sample inside of a space,
#' calculated prior from a compendium of known samples (Check ?SpaceMaker for more info).
#'
#' @param InputData A matrix of input gene expressions to be analyzed, with genes as rows and samples as columns. Corresponding Entrez Gene
#' IDs must be assigned to 'rownames' of the matrix, and name of each sample/column should be written in 'colnames'.
#' @param Space
#' @return NULL
#'
#' @examples print()
#'
#' @export calculatePhysioMap
#Pre-function for dispaching the calculating to the write function:
calculatePhysioMap <- function(InputData, Space, GenesRatio = 0.05,
                               PARALLEL = FALSE, TTEST = FALSE,
                               WilxSTATICResponse = F, NumbrOfCores=NA){
  UseMethod("calculatePhysioMap", object = Space)
}

#' @export calculatePhysioMap.default
##The main function for calculating PhysioScores:
calculatePhysioMap.default <- function(InputData, Space, GenesRatio = 0.05,
                                       PARALLEL = FALSE, TTEST = FALSE,
                                       WilxSTATICResponse = F, NumbrOfCores=NA){
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
    foreach(i = 1:NSamples, .combine='cbind', .final=as.matrix, .export=c("TTestWrapper","wilTestWrapper"), .packages = "PhysioSpaceMethods") %dopar% {
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
      apply(X = Space, MARGIN = 2, FUN = if(TTEST) TTestWrapper else wilTestWrapper,
            iplus=iplus, iminus = iminus, WilxSTATICResponse = WilxSTATICResponse) # this apply can be written as simple for loop also, results
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
                                    PARALLEL = FALSE, TTEST = FALSE,
                                    WilxSTATICResponse = F, NumbrOfCores=NA){

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
