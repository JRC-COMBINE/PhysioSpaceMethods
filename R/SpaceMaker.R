#' @title Creates PhysioSpaces
#'
#' @description This function uses 'Big Data' to make robust
#' 'Physiological Vectors' in N dimensional spaces, within which you
#' can map new data to extract physiological information from a new data set.
#'
#' @param GeneExMatrix A matrix of input gene expressions or a
#' SummarizedExperiment object, based on which the Physiological
#' Space is made.
#'
#' In case of a matrix, GeneExMatrix is supposed to have genes
#' as rows and samples as columns. Corresponding Entrez Gene
#' IDs must be assigned to 'rownames', and name of each sample
#' should be written in 'colnames' of the matrix.
#'
#' In case of a SummarizedExperiment object, GeneExMatrix
#' must have a component named 'EntrezID' in its rowData.
#' It is also expected (but not mandatory) for GeneExMatrix
#' to have a component named 'SampleName' in its colData.
#' The gene expressions in GeneExMatrix is extracted by the
#' function assay(), meaning in case GeneExMatrix contains
#' multiple assays, only the first one is used.
#'
#' Unless 'DESIGN' and 'CONTRASTs' inputs are provided by the
#' user, spaceMaker supposes the label of the first column
#' (colnames(GeneExMatrix)[1]) to be the reference of the
#' experiment and uses all the samples with this label as control.
#'
#' @param DESIGN (Optional) Design matrix of GeneExMatrix, made by
#' the function model.matrix(). If it's not provided, spaceMaker()
#' will make a design matrix based on sample names of GeneExMatrix.
#'
#' @param CONTRASTs (Optional) character vector or list specifying
#' contrasts. If it's not provided, spaceMaker()
#' will make the CONTRASTs with the assumption that
#' sample names of first column is the label of the control or reference.
#' REMEMBER that expected user-defined CONTRASTs format changes based
#' on the LinearOrRNASeq input: in case LinearOrRNASeq='Linear',
#' CONTRASTs is expected to work as an input for
#' limma::makeContrasts(). And when LinearOrRNASeq='RNASeq',
#' CONTRASTs is used as an input for DESeq2::results().
#'
#' @param Output A character specifying the output format of
#' spaceMaker(). The default value is 'PhysioScore', which will
#' return -log2(p value)*sign(fold change). It is also possible
#' to obtain fold change by Output='FoldChange', or obtain
#' the fitted model by having Output ='Model'.
#'
#' @param LinearOrRNASeq A character which determines what type of
#' modelling is ought to be used when making the PhysioSpace.
#' If it's possible to do linear modelling on the data, e.g. data
#' is log normal micro-array gene expression data or
#' limma::voom-transformed RNA-seq data, then LinearOrRNASeq should
#' be 'Linear'. In this case limma package is used in the
#' calculations. But in case your GeneExMatrix input is an RNA-seq
#' count matrix, you should pass 'RNASeq' to LinearOrRNASeq. In
#' this case DESeq2 package is used for calculations.
#'
#'
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#'
#' @return Depending on the 'Output' argument, the returned value is
#' either a matrix, or a model. If Output = "PhysioScore",
#' a matrix is returned, with genes in rows and Physiological axes
#' on the columns. In this case, values inside this matrix are
#' PhysioScores (-log2(p value)*sign(fold change)). In case of
#' Output = "FoldChange", a matrix of fold changes is returned. And
#' if Output = "Model", the fitted model by limma::lmFit() or
#' DESeq2::DESeq() is returned. REMEMBER that when user provides
#' 'DESIGN' input argument, colnames of the returned matrix remains
#' empty and are needed to be assigned by the user.
#'
#'
#' @examples
#'  INPTMat <-
#'    matrix(
#'      data = rnorm(n = 180000, mean = 8, sd = 6),
#'      nrow = 20000,
#'      dimnames = list(paste0("g", 1:20000), c(
#'        rep("Ctrl", 3), rep("Cancer1", 3), rep("Cancer2", 3)
#'      ))
#'    ) #Simulated DNA-array gene expression matrix
#'  LinearSpaceOfINPTMat <- spaceMaker(GeneExMatrix = INPTMat)
#'
#'  INPTMatRNASeq <-
#'    matrix(
#'      data = rnbinom(n = 180000, size = 1.5, prob = 0.01),
#'      nrow = 20000,
#'      dimnames = list(paste0("g", 1:20000), c(
#'        rep("Ctrl", 3), rep("Cancer1", 3), rep("Cancer2", 3)
#'      ))
#'    ) #Simulated RNA-seq gene expression matrix
#'  NotLinearSpaceOfINPTMatRNASeq <-
#'    spaceMaker(GeneExMatrix = INPTMatRNASeq, LinearOrRNASeq = "RNASeq")
#'
#'
#'  library(SummarizedExperiment)
#'  INPTMat_SE <- SummarizedExperiment(
#'      assays = list(GEX = INPTMat),
#'      rowData = data.frame("EntrezID" = rownames(INPTMat)),
#'      colData = data.frame("SampleName" = colnames(INPTMat))
#'  ) #Simulated DNA-array gene expression SummarizedExperiment obj.
#'  LinearSpaceOfINPTMat_SE <- spaceMaker(GeneExMatrix = INPTMat_SE)
#'
#'
#' @export
spaceMaker <- function(GeneExMatrix, DESIGN = NA, CONTRASTs = NA,
                       Output = "PhysioScore", LinearOrRNASeq = "Linear"){

    if(!any(LinearOrRNASeq == c("Linear","RNASeq"))){
        stop("'LinearOrRNASeq' input should be either 'Linear' or 'RNASeq'")
    }
    UserDefinedDesign <- TRUE
    GeneExMatrix <- .inptPreparer(InputData = GeneExMatrix)

    # 'Design' construction from colnames of GeneExMatrix:
    if(is.na(DESIGN)){
        DESIGN <- model.matrix(object =
                                   ~ 0 + factor(x = colnames(GeneExMatrix),
                                                levels = unique(colnames(GeneExMatrix))))
        #made a factor with first column as first level
        colnames(DESIGN) <- unique(colnames(GeneExMatrix))

        UserDefinedDesign <- FALSE
    }

    # 'Contrast' construction from colnames of GeneExMatrix:
    if(LinearOrRNASeq=="Linear"){
        if(is.na(CONTRASTs)){
            #Making the 'contrast's from colnames of GeneExMatrix:
            CONTRASTs <- paste(colnames(DESIGN)[-1],
                               colnames(DESIGN)[1], sep = "-")
        }
    } else {
        #Making colData for DESeqDataSetFromMatrix:
        colDataForDESeqModel = data.frame("CONDITION" = colnames(GeneExMatrix))
        colnames(GeneExMatrix) <- make.names(colnames(GeneExMatrix),
                                             unique = TRUE)
        ##Because DESeqDataSetFromMatrix forces colnames of countData as
        #rownames to colData, and colData is a data.frame so repeat
        #in those names mean error.
        if(is.na(CONTRASTs)){
            #Making the 'contrast's from colnames of GeneExMatrix:
            CONTRASTs <- list()
            for(K in 2:ncol(DESIGN)){
                CONTRASTs[[colnames(DESIGN)[K]]] <-
                    list(colnames(DESIGN)[K],
                         colnames(DESIGN)[1])
            }
        }
    }

    #MAIN model fitting:
    if(LinearOrRNASeq=="Linear"){
        MODELStructured <- lmFit(GeneExMatrix, DESIGN)
        cont.matrix <- makeContrasts(contrasts = CONTRASTs, levels=DESIGN)
        MODELCalculated <- contrasts.fit(MODELStructured, cont.matrix)
        MODELCalculated <- eBayes(MODELCalculated)
    } else {
        MODELStructured <- DESeqDataSetFromMatrix(countData = GeneExMatrix,
                                                  colData = colDataForDESeqModel,
                                                  design = DESIGN)
        MODELCalculated <- DESeq(MODELStructured, quiet = TRUE)
    }

    #Extracting aimed results:
    MODELResults <- list()
    pb <- progress_bar$new(format = "(:spin) [:bar] :percent eta: :eta",
                           total = length(CONTRASTs), clear = FALSE)
    for(Const in CONTRASTs){
        MODELResults[[length(MODELResults) + 1]] <-
            if(LinearOrRNASeq=="Linear"){
                topTable(fit = MODELCalculated, coef = Const,
                         adjust.method="BH",
                         number = Inf, sort.by = "none")
            } else {
                results(object = MODELCalculated,
                        pAdjustMethod = "BH",
                        contrast = Const)
            }
        pb$tick()
    }

    #(Calculating and) Returning the preferred variables:
    PVALs <- if(LinearOrRNASeq=="Linear") "adj.P.Val" else "padj"
    LFCs <- if(LinearOrRNASeq=="Linear") "logFC" else "log2FoldChange"
    if(Output == "PhysioScore"){
        Outi <- vapply(X = MODELResults,
                       function(x) -log2(x[[PVALs]])*sign(x[[LFCs]]),
                       FUN.VALUE = numeric(length = nrow(GeneExMatrix)))
    } else if(Output == "FoldChange"){
        Outi <- vapply(X = MODELResults, function(x) x[[LFCs]],
                       FUN.VALUE = numeric(length = nrow(GeneExMatrix)))
    } else if(Output == "Model"){
        return(MODELStructured)
    } else {
        stop("'Output' value is incorrect.")
    }
    rownames(Outi) <- rownames(GeneExMatrix)
    if(!UserDefinedDesign) colnames(Outi) <- colnames(DESIGN)[-1]
    return(Outi)

}

