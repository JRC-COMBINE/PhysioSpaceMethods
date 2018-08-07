#' @title Mapping new data into a physiological-space
#'
#' @description calculatePhysioMap computes mapped values of each input sample
#' inside of a space,
#' calculated prior from a compendium of known samples.
#'
#' @param InputData A matrix of input gene expressions to be analyzed, with
#' genes as rows and samples as columns. Corresponding Entrez Gene
#' IDs must be assigned to 'rownames' of the matrix, and name of each
#' sample/column should be written in 'colnames'. REMEMBER that the gene
#' expressions in 'InputData' should be relative; e.g. fold change or signed
#' p value of a statistical test.
#' In case user has their own list of significantly up and down regulated
#' genes, it is also possible for InputData to be a list, containing
#' Entrez IDs (or any other identifier which is used as rownames
#' in 'Space') of up regulated genes in InputData[[1]] and
#' Entrez IDs (or any other identifier which is used as rownames
#' in 'Space') of down regulated genes in InputData[[2]].
#' Having a list InputData is usually slow and restrictive, hence,
#' mainly it is not recommended.
#' @param Space The space in which the 'InputData' will be mapped. Just as
#' 'InputData', it should be a matrix with genes as rows and samples as
#' columns, with corresponding Entrez Gene IDs in 'rownames' of the matrix,
#' and name of each axis of the space written in 'colnames'.
#' @param GenesRatio The ratio of gene expression values to be considered in
#' the calculation. In high dimensional omics data, signal to noise ratio
#' has a direct relation with the relative magnitude of expressions.
#' We aim to remove the noisy genes, hence we only keep the "GenesRatio*100"
#' percent highest and lowest gene expression values of each sample.
#' GenesRatio should be a numerical value between 0 and 1. Default value
#' is 0.05. This input is only applicable when InputData is a matrix.
#' @param PARALLEL Logical value indicating if calculation should be done in
#' parallel. Default value is FALSE.
#' It is not recommended to use PARALLEL=TRUE on small datasets since due to
#' large overhead, it could take more time than using PARALLEL=FALSE.
#' This input is only applicable when InputData is a matrix.
#' @param NumbrOfCores Number of cores to be used when 'PARALLEL' is TRUE.
#' Default is NA which will result in the program using
#' 'all available cores - 1'. Assigning a number higher than
#' parallel::detectCores() will result in an error.
#' @param TTEST Logical value indicating if t.test should be done in place
#' of the default wilcoxon rank-sum test (more info can be found
#' in the original PhysioSpace: Lenz et. al., PLOS One 2013). Using t.test
#' will speed up calculations. Default value is FALSE.
#' @param STATICResponse Logical value indicating if 'statistic' should be
#' returned rather than the default 'signed p value'. Default value is FALSE.
#' @param ImputationMethod Imputation method to use in case of missing
#' values. Available methods are "PCA" and "KNN". Default is "PCA".
#' This input is only applicable when InputData is a matrix.
#' @param ParallelMethod Parallel method to use when PARALLEL=TRUE.
#' Two methods are implemented so far: "parCapply"
#' which uses parCapply function of parallel package, and "foreach" which
#' uses foreach package to process in parallel.
#' Speed-wise, foreach has an edge on parCapply, but the latter
#' is more stable. Hence, we recommend to use parCapply.
#' The default value for ParallelMethod is "parCapply".
#' This input is only applicable when InputData is a matrix.
#'
#' @import progress
#'
#' @details PhysioSpace is a robust
#' statistical method for relating high dimensional omics data sets from
#' heterogeneous sources using shared physiological
#' processes. It is designed to take advantage of the vast availability of
#' public omics data, which in combination with
#' statistical approaches makes a potent tool capable of analyzing
#' heterogeneous biological data sets.
#' 'calculatePhysioMap' is the main analytical function of the package.
#' It uses a nonlinear mapping function to relate the
#' unknown input data with a physiological space. Physiological spaces
#' are mathematical spaces build upon known
#' physiological data, using the 'spaceMaker' function.
#'
#' @return Matrix of mapped 'InputData' values in 'Space', with rows
#' corresponding to axes of 'Space' and columns representing
#' samples in 'InputData'. Mapped values are signed p value when
#' STATICResponse==FALSE, and are 'statistic' value when
#' STATICResponse==TRUE (more info can be found in the original
#' PhysioSpace paper: Lenz et. al., PLOS One 2013).
#'
#' @references Lenz, M., Schuldt, B. M., Müller, F. J., & Schuppert,
#' A. (2013). PhysioSpace: relating gene expression experiments from
#' heterogeneous sources using shared physiological processes.
#' PLoS One, 8(10), e77627.
#'
#' @examples
#'  SimulatedGeneExpressionData <- matrix(
#'    rnorm(n = 100000, mean = 0,
#'          sd = 100),
#'    ncol = 10,
#'    dimnames = list(1:10000, 1:10)
#'  )
#'  SimulatedReferenceSpace <- matrix(
#'    rnorm(n = 100000, mean = 0,
#'          sd = 100),
#'    ncol = 10,
#'    dimnames = list(1:10000, 11:20)
#'  )
#'  calculatePhysioMap(InputData = SimulatedGeneExpressionData,
#'                     Space = SimulatedReferenceSpace)
#'  if (parallel::detectCores() > 1) {
#'    #More than one core is needed for parallel processing
#'    calculatePhysioMap(
#'      InputData = SimulatedGeneExpressionData,
#'      Space = SimulatedReferenceSpace,
#'      PARALLEL = TRUE,
#'      NumbrOfCores = 2,
#'      GenesRatio = 0.01,
#'      STATICResponse = FALSE,
#'      TTEST = TRUE
#'    )
#'  }
#'
#' @export
#Pre-function for dispaching the calculating to the write function:
calculatePhysioMap <- function(InputData, Space, GenesRatio = 0.05,
                                PARALLEL = FALSE, NumbrOfCores=NA,
                                TTEST = FALSE, STATICResponse = FALSE,
                                ImputationMethod = "PCA",
                                ParallelMethod = "parCapply"){
    UseMethod("calculatePhysioMap")
}


##The main function for calculating PhysioScores:
#' @export
calculatePhysioMap.default <- function(InputData, Space, GenesRatio = 0.05,
                                        PARALLEL = FALSE, NumbrOfCores=NA,
                                        TTEST = FALSE, STATICResponse = FALSE,
                                        ImputationMethod = "PCA",
                                        ParallelMethod = "parCapply"){
    #Initializing:
    inputChecker(InputData, Space, ImputationMethod)
    NGenes <- nrow(InputData)
    NSamples <- ncol(InputData)
    physioMap <- matrix(NA, ncol(Space), NSamples)
    if(PARALLEL) cl <- parallelInitializer(NumbrOfCores=NumbrOfCores)
    pb <- progress_bar$new(format = "(:spin) [:bar] :percent eta: :eta",
                            total = ifelse(PARALLEL,NSamples/length(cl),
                                            NSamples),
                            clear = FALSE)
    #
    #Main:
    if(PARALLEL) {
        if(ParallelMethod == "parCapply"){
            physioMap <- calculatePhysioMapCore_parCapply(InputData,
                                                            Space,
                                                            NSamples,
                                                            GenesRatio,
                                                            NGenes,
                                                            STATICResponse,
                                                            pb,
                                                            TTEST,
                                                            cl)
        } else if(ParallelMethod == "foreach"){
            physioMap <- calculatePhysioMapCore_foreach(InputData, Space,
                                                        NSamples, GenesRatio,
                                                        NGenes, STATICResponse,
                                                        pb, TTEST)
        } else {
            stop("The specified ParallelMethod is not implemented!")
        }
        stopCluster(cl)
    } else {
        physioMap <- apply(X = InputData,
                            MARGIN = 2, FUN = singleThreadOfPhysioCalc,
                            Space=Space, GenesRatio=GenesRatio,
                            NGenes=NGenes, STATICResponse=STATICResponse,
                            pb=pb, TTEST=TTEST)
    }

    rownames(physioMap) = colnames(Space)
    colnames(physioMap) = colnames(InputData)
    return(physioMap)
}


##PhysioSpace for gene lists as a Input:
#' @export
calculatePhysioMap.list <- function(InputData, Space, GenesRatio = 0.05,
                                    PARALLEL = FALSE, NumbrOfCores=NA,
                                    TTEST = FALSE, STATICResponse = FALSE,
                                    ImputationMethod = "PCA",
                                    ParallelMethod = "parCapply"){
    UpIndx <- na.omit(match(InputData[[1]],rownames(Space)))
    DownIndx <- na.omit(match(InputData[[2]],rownames(Space)))
    if(length(UpIndx) < 2){
        stop(paste("Not enough up regulated genes were",
                    "found in rownames of Space"))
    }
    if(length(DownIndx) < 2) {
        stop(paste("Not enough down regulated genes",
                    "were found in rownames of Space"))
    }
    as.matrix(apply(X = Space, MARGIN = 2,
                    FUN = if(TTEST) tTestWrapper else wilTestWrapper,
                    iplus=UpIndx, iminus = DownIndx,
                    STATICResponse = STATICResponse))
}
