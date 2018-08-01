context("Test basic Physio Score calculations")

SimulatedGeneExpressionData <- matrix(rnorm(n = 100000, mean = 0, sd = 100),ncol = 10, dimnames = list(1:10000,1:10))
SimulatedReferenceSpace <- matrix(rnorm(n = 100000, mean = 0, sd = 100),ncol = 10, dimnames = list(1:10000,11:20))
SelfSCORES <- calculatePhysioMap(InputData = SimulatedGeneExpressionData, SimulatedGeneExpressionData)
SCORES <- calculatePhysioMap(InputData = SimulatedGeneExpressionData, SimulatedReferenceSpace)

expectedSelfMax <- -log2(wilcox.test(x = 1:round(nrow(SimulatedGeneExpressionData)*0.05),
                               y = (round(nrow(SimulatedGeneExpressionData)*0.05)+1):
                                 (2*round(nrow(SimulatedGeneExpressionData)*0.05)))$p.value)

test_that("'calculatePhysioMap' has to have a matrices as input and space",{
  expect_is(SelfSCORES,"matrix")
  expect_length(SelfSCORES, ncol(SimulatedGeneExpressionData)*ncol(SimulatedReferenceSpace))
  expect_equal(max(SelfSCORES),expectedSelfMax)
  expect_lt(max(SCORES),expectedSelfMax)
})
