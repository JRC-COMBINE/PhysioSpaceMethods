## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = ">"
)

## ----message=FALSE------------------------------------------------------------
library(SummarizedExperiment) #SummarizedExperiment is needed for 
                            #working with RangedSummarizedExperiment objects.
library(EnrichmentBrowser)  # For converting IDs of a
                            # RangedSummarizedExperiment object.
library(org.Hs.eg.db) #For ID conversion
library(ExperimentHub) #Downloads datasets from ExperimentHub web service
library(PhysioSpaceMethods) #Main package

## ----message=FALSE------------------------------------------------------------
hub <- ExperimentHub()
HumanAffys <- query(hub, "HumanAffyData")
LukkData <- HumanAffys[["EH177"]]

## -----------------------------------------------------------------------------
LukkDataMatrix <- exprs(LukkData)

## -----------------------------------------------------------------------------
colnames(LukkDataMatrix) <- make.names(pData(LukkData)$Groups_369)
#'make.names' is used so that some specific names won't break the pipeline
#later on. 

## -----------------------------------------------------------------------------
LukkDataMatrix <- cbind(apply(LukkDataMatrix,1,mean),LukkDataMatrix)
colnames(LukkDataMatrix)[1] <- "Ctrl"

## ----message=FALSE------------------------------------------------------------
LukkSpace <- spaceMaker(GeneExMatrix = LukkDataMatrix)

## -----------------------------------------------------------------------------
#Download:
load(url(paste0("https://www.ebi.ac.uk/gxa/experiments-content/",
                "E-MTAB-2836/static/E-MTAB-2836-atlasExperimentSummary.Rdata")))

## -----------------------------------------------------------------------------
#Making the gene expression matrix:
EMTAB2836CountMatrix <- assay(experimentSummary$rnaseq)

## ----message=FALSE------------------------------------------------------------
#Converting Ensembl to Entrez IDs:
ENSEMBL2EG <- as.list(org.Hs.egENSEMBL2EG)
IDIndx <- match(rownames(EMTAB2836CountMatrix), names(ENSEMBL2EG), nomatch = 0)
EMTAB2836CountMatrix <- EMTAB2836CountMatrix[IDIndx!=0,]
rownames(EMTAB2836CountMatrix) <- sapply(ENSEMBL2EG[IDIndx], function(x) x[1])

## -----------------------------------------------------------------------------
#Assigning colnames:
colnames(EMTAB2836CountMatrix) <-
                            colData(experimentSummary$rnaseq)$organism_part

## -----------------------------------------------------------------------------
#Calculating Fold-Changes:
EMTAB2836CountMatrixRelativ <- EMTAB2836CountMatrix -
                                        apply(EMTAB2836CountMatrix,1,mean)

## -----------------------------------------------------------------------------
#Loading the object into R:
EMTAB2836_SEObj <- experimentSummary$rnaseq
#
#Converting Ensembl to Entrez IDs:
EMTAB2836_SEObj <- idMap(EMTAB2836_SEObj, org="hsa",
                                                from="ENSEMBL", to="ENTREZID")
#Making EntrezID in rowData:
rowData(EMTAB2836_SEObj) <- data.frame("EntrezID" =
                                           rownames(EMTAB2836_SEObj))

## ----message=FALSE------------------------------------------------------------
EMTAB2836_SEObj <- experimentSummary$rnaseq
#
#Converting Ensembl to Entrez IDs, and making EntrezID in rowData:
ENSEMBL2EG <- as.list(org.Hs.egENSEMBL2EG)
IDIndx <- match(rownames(EMTAB2836_SEObj), names(ENSEMBL2EG), nomatch = 0)
EMTAB2836_SEObj <- EMTAB2836_SEObj[IDIndx!=0,]
rowData(EMTAB2836_SEObj) <- data.frame("EntrezID" =
              sapply(ENSEMBL2EG[IDIndx], function(x) x[1]))

## -----------------------------------------------------------------------------
#Assigning SampleName:
names(colData(EMTAB2836_SEObj))[names(colData(EMTAB2836_SEObj)) ==
                                    "organism_part"] <- "SampleName"

## -----------------------------------------------------------------------------
#Calculating Fold-Changes:
assay(EMTAB2836_SEObj) <- assay(EMTAB2836_SEObj) -
                                        apply(assay(EMTAB2836_SEObj),1,mean)

## ----message=FALSE------------------------------------------------------------
#Choosing 5 random samples:
set.seed(seed = 0) #So results would be reproducible
Samples5Random <- sample(x = 1:ncol(EMTAB2836CountMatrixRelativ), size = 5)
#Main calculation:
RESULTS5 <- calculatePhysioMap(
    InputData = EMTAB2836CountMatrixRelativ[,Samples5Random],
    Space = LukkSpace)

## ----message=FALSE------------------------------------------------------------
#Main calculation, for SummarizedExperiment input obj.:
RESULTS5_SE <- calculatePhysioMap(
    InputData = EMTAB2836_SEObj[,Samples5Random],
    Space = LukkSpace)

## -----------------------------------------------------------------------------
#Check to see if the results are the same, regardless of input format:
identical(RESULTS5, RESULTS5_SE)

## ----message=FALSE------------------------------------------------------------
#Main calculation in parallel:
RESULTS5 <- calculatePhysioMap(
    InputData = EMTAB2836CountMatrixRelativ[,Samples5Random],
    Space = LukkSpace, NumbrOfCores = 2)

## ----fig.width=10, fig.asp=1--------------------------------------------------
#Plotting the results:
PhysioHeatmap(PhysioResults = RESULTS5, main = "RNA-seq vs Microarray",
            SymmetricColoring = TRUE, SpaceClustering = FALSE,
            Space = LukkSpace, ReducedPlotting = 5)

## -----------------------------------------------------------------------------
sessionInfo()

