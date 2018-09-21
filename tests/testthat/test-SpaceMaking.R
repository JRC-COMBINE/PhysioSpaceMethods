context("Test basic Space making")

INPTMat <-
    matrix(
        data = rnorm(n = 180000, mean = 8, sd = 6),
        nrow = 20000,
        dimnames = list(paste0("g", 1:20000), c(
            rep("Ctrl", 3), rep("Cancer1", 3), rep("Cancer2", 3)
        ))
    ) #Simulated DNA-array gene expression matrix
LinearSpaceOfINPTMat <- spaceMaker(GeneExMatrix = INPTMat,
                                    Output = "PhysioScore")
LinearSpaceOfINPTMatFD <- spaceMaker(GeneExMatrix = INPTMat,
                                    Output = "FoldChange")
LinearSpaceOfINPTMatModl <- spaceMaker(GeneExMatrix = INPTMat,
                                    Output = "Model")

INPTMatRNASeq <-
    matrix(
        data = rnbinom(n = 180000, size = 1.5, prob = 0.01),
        nrow = 20000,
        dimnames = list(paste0("g", 1:20000), c(
            rep("Ctrl", 3), rep("Cancer1", 3), rep("Cancer2", 3)
        ))
    ) #Simulated RNA-seq gene expression matrix
NotLinearSpaceOfINPTMatRNASeq <-
    spaceMaker(GeneExMatrix = INPTMatRNASeq,
                LinearOrRNASeq = "RNASeq", Output = "PhysioScore")
NotLinearSpaceOfINPTMatRNASeqFD <-
    spaceMaker(GeneExMatrix = INPTMatRNASeq,
                LinearOrRNASeq = "RNASeq", Output = "FoldChange")
NotLinearSpaceOfINPTMatRNASeqModl <-
    spaceMaker(GeneExMatrix = INPTMatRNASeq,
                LinearOrRNASeq = "RNASeq", Output = "Model")

library(SummarizedExperiment)
INPTMat_SE <- SummarizedExperiment(
    assays = list(GEX = INPTMat),
    rowData = data.frame("EntrezID" = rownames(INPTMat)),
    colData = data.frame("SampleName" = colnames(INPTMat))
) #Simulated DNA-array gene expression SummarizedExperiment obj.
LinearSpaceOfINPTMat_SE <- spaceMaker(GeneExMatrix = INPTMat_SE)


test_that("'spaceMaker' returns a specific matrix",{
    expect_is(LinearSpaceOfINPTMat,"matrix")
    expect_equal(nrow(LinearSpaceOfINPTMat), nrow(INPTMat))
    expect_is(LinearSpaceOfINPTMatFD,"matrix")
    expect_equal(nrow(LinearSpaceOfINPTMatFD), nrow(INPTMat))
    expect_is(NotLinearSpaceOfINPTMatRNASeq,"matrix")
    expect_equal(nrow(NotLinearSpaceOfINPTMatRNASeq), nrow(INPTMat))
    expect_is(NotLinearSpaceOfINPTMatRNASeqFD,"matrix")
    expect_equal(nrow(NotLinearSpaceOfINPTMatRNASeqFD), nrow(INPTMat))
    expect_is(LinearSpaceOfINPTMat_SE,"matrix")
    expect_equal(nrow(LinearSpaceOfINPTMat_SE), nrow(INPTMat))
})

test_that("'spaceMaker' could return a model",{
    expect_is(LinearSpaceOfINPTMatModl, "MArrayLM")
    expect_is(NotLinearSpaceOfINPTMatRNASeqModl, "DESeqDataSet")
})

test_that("'spaceMaker' scores are bounded",{
    expect_lte(max(2^-abs(LinearSpaceOfINPTMat)), 1)
    expect_gte(min(2^-abs(LinearSpaceOfINPTMat)), 0)
    expect_lte(max(2^-abs(NotLinearSpaceOfINPTMatRNASeq),
                   na.rm = TRUE), 1)
    expect_gte(min(2^-abs(NotLinearSpaceOfINPTMatRNASeq),
                   na.rm = TRUE), 0)
    expect_lte(max(2^-abs(LinearSpaceOfINPTMat_SE)), 1)
    expect_gte(min(2^-abs(LinearSpaceOfINPTMat_SE)), 0)
})
