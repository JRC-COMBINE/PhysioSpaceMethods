context("Test gene expression imputation")

MatToImpute <-
    matrix(
        rnorm(n = 10000, mean = 0, sd = 100),
        ncol = 10,
        dimnames = list(1:1000, 1:10)
    )
MatToImpute[sample(x = 1:length(MatToImpute),
                   size = length(MatToImpute) / 20)] <- NA
ImputedMatPCA <-
    .imputeMissingGeneExpression(InptGEX = MatToImpute, METHOD = "PCA")
ImputedMatKNN <-
    .imputeMissingGeneExpression(InptGEX = MatToImpute, METHOD = "KNN")

test_that("imputation returns a non-na matrix",{
    expect_is(ImputedMatPCA,"matrix")
    expect_length(ImputedMatPCA, ncol(MatToImpute)*nrow(MatToImpute))
    expect_false(anyNA(ImputedMatPCA))
    expect_is(ImputedMatKNN,"matrix")
    expect_length(ImputedMatKNN, ncol(MatToImpute)*nrow(MatToImpute))
    expect_false(anyNA(ImputedMatKNN))
})

test_that("Unknown method should make an error",{
    expect_error(.imputeMissingGeneExpression(InptGEX =
                                                 MatToImpute,
                                             METHOD = "NewMethod"))
})
