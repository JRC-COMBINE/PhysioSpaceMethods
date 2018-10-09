context("Test input preparation")

SimulatedGeneExpressionData <- matrix(
    rnorm(n = 100000, mean = 0,
          sd = 100),
    ncol = 10,
    dimnames = list(1:10000, 1:10)
)
SimulatedGeneExpressionData_checked <-
    .inptPreparer(SimulatedGeneExpressionData)

library(SummarizedExperiment)
SimulatedGeneExpressionData_SE <- SummarizedExperiment(
    assays = list(GEX = SimulatedGeneExpressionData),
    rowData = data.frame("EntrezID" =
                             rownames(SimulatedGeneExpressionData)),
    colData = data.frame("SampleName" =
                             colnames(SimulatedGeneExpressionData))
)
SimulatedGeneExpressionData_SE_checked <-
    .inptPreparer(SimulatedGeneExpressionData_SE)

test_that("preparation returns a matrix",{
    expect_is(SimulatedGeneExpressionData_checked,"matrix")
    expect_is(SimulatedGeneExpressionData_SE_checked,"matrix")
    expect_length(SimulatedGeneExpressionData_checked, 10*10000)
    expect_length(SimulatedGeneExpressionData_SE_checked, 10*10000)
    expect_equal(SimulatedGeneExpressionData_checked,
                 SimulatedGeneExpressionData)
    expect_equal(SimulatedGeneExpressionData_checked,
                 SimulatedGeneExpressionData_SE_checked)
})

test_that("Unknown object should make an error",{
    expect_error(.inptPreparer(InputData = "Some random Character string"))
})
