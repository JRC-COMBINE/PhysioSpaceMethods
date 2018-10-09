context("Test inptChecking")

SimulatedGeneExpressionData <-
    matrix(
        rnorm(n = 100000, mean = 0, sd = 100),
        ncol = 10,
        dimnames = list(1:10000, 1:10)
    )
SimulatedGeneExpressionData[sample(x = 1:length(SimulatedGeneExpressionData),
                    size = length(SimulatedGeneExpressionData) / 20)] <- NA
.inptChecker(InputData = SimulatedGeneExpressionData[, 1:5],
            Space = SimulatedGeneExpressionData[sample(1:10000), 6:10])


test_that("'.inptChecker' returns a 'InputData' to check",{
    expect_is(InputData,"matrix")
    expect_length(InputData, 5*10000)
})

test_that("'.inptChecker' also returns a 'Space' to check",{
    expect_is(Space,"matrix")
    expect_length(Space, 5*10000)
    expect_equal(dim(InputData),dim(Space))
})
