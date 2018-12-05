context("Test basic Physio Score calculations")

SimulatedGeneExpressionData <- matrix(rnorm(n = 10000, mean = 0, sd = 100),
                                    ncol = 10, dimnames = list(1:1000,1:10))
SimulatedReferenceSpace <- matrix(rnorm(n = 10000, mean = 0, sd = 100),
                                    ncol = 10, dimnames = list(1:1000,11:20))
SelfSCORES <- calculatePhysioMap(InputData = SimulatedGeneExpressionData,
                                    SimulatedGeneExpressionData)
SCORES <- calculatePhysioMap(InputData = SimulatedGeneExpressionData,
                                SimulatedReferenceSpace)
SCORESParallel <- calculatePhysioMap(InputData = SimulatedGeneExpressionData,
                                        Space = SimulatedReferenceSpace,
                                        NumbrOfCores = 2)
SCORESParallel2 <- calculatePhysioMap(InputData = SimulatedGeneExpressionData,
                                        Space = SimulatedReferenceSpace,
                                        NumbrOfCores = MulticoreParam(2))
SCORESParallel3 <- calculatePhysioMap(InputData = SimulatedGeneExpressionData,
                                        Space = SimulatedReferenceSpace,
                                        NumbrOfCores = SnowParam(2))

expectedSelfMax <-
    -log2(wilcox.test(
        x = 1:round(nrow(SimulatedGeneExpressionData) * 0.05),
        y = (round(nrow(
            SimulatedGeneExpressionData
        ) * 0.05) + 1):(2 * round(nrow(
            SimulatedGeneExpressionData
        ) * 0.05))
    )$p.value)

test_that("'calculatePhysioMap' has to have a matrices as input and space",{
  expect_is(SelfSCORES,"matrix")
  expect_length(SelfSCORES,
                ncol(SimulatedGeneExpressionData)*ncol(SimulatedReferenceSpace))
  expect_equal(max(SelfSCORES),expectedSelfMax)
  expect_lt(max(SCORES),expectedSelfMax)
  expect_error(calculatePhysioMap(InputData = SimulatedGeneExpressionData,
                                  Space = SimulatedReferenceSpace,
                                  NumbrOfCores = 2, ParallelMethod = "NewMeth"))
  expect_equal(SCORES,SCORESParallel)
  expect_equal(SCORES,SCORESParallel2)
  expect_equal(SCORES,SCORESParallel3)
})

test_that("'calculatePhysioMap' should work with single input or space",{

    SSSingle <- calculatePhysioMap(InputData =
                                SimulatedGeneExpressionData[,1,drop=FALSE],
                                     SimulatedGeneExpressionData)
    SSingle <- calculatePhysioMap(InputData = SimulatedGeneExpressionData,
                                 SimulatedReferenceSpace[,1,drop=FALSE])
    SSPSingle <- calculatePhysioMap(InputData =
                                    SimulatedGeneExpressionData[,1,drop=FALSE],
                            Space = SimulatedReferenceSpace[,1,drop=FALSE],
                            NumbrOfCores = 2)
    expect_is(SSSingle,"matrix")
    expect_is(SSingle,"matrix")
    expect_length(SSSingle, ncol(SimulatedReferenceSpace))
    expect_equal(SSingle[1],as.numeric(SSPSingle))
})
