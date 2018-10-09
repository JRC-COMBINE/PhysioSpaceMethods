context("Test inptImputing")

SimulatedGeneExpressionData <-
    matrix(
        rnorm(n = 100000, mean = 0, sd = 100),
        ncol = 10,
        dimnames = list(1:10000, 1:10)
    )
SimulatedGeneExpressionData[sample(x = 1:length(SimulatedGeneExpressionData),
                    size = length(SimulatedGeneExpressionData) / 20)] <- NA
SimulatedGeneExpressionDataImputed <-
    .imputeMissingGeneExpression(SimulatedGeneExpressionData,
                                 METHOD="PCA")
SimulatedGeneExpressionDataImputedKNN <-
    .imputeMissingGeneExpression(SimulatedGeneExpressionData,
                                 METHOD="KNN")


test_that("'.imputeMissingGeneExpression' returns a matrix with imputed data",{
    expect_is(SimulatedGeneExpressionDataImputed,"matrix")
    expect_false(anyNA(SimulatedGeneExpressionDataImputed))
    expect_equal(dim(SimulatedGeneExpressionDataImputed),
                 dim(SimulatedGeneExpressionData))

    expect_is(SimulatedGeneExpressionDataImputedKNN,"matrix")
    expect_false(anyNA(SimulatedGeneExpressionDataImputedKNN))
    expect_equal(dim(SimulatedGeneExpressionDataImputedKNN),
                 dim(SimulatedGeneExpressionData))
})
