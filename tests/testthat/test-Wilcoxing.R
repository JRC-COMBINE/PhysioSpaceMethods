context("Test Wilcoxon rank-sum wrapper function used by calculatePhysioMap")

SimulatedReferenceSpace <- matrix(rnorm(n = 100000, mean = 0, sd = 100),
                                  ncol = 10, dimnames = list(1:10000,11:20))

loggedPs <- .wilTest(
        ReferencesJ = SimulatedReferenceSpace[, 4],
        iplus = sample(
          1:nrow(SimulatedReferenceSpace),
          size = nrow(SimulatedReferenceSpace) / 20
        ),
        iminus = sample(
          1:nrow(SimulatedReferenceSpace),
          size = nrow(SimulatedReferenceSpace) / 20
        ),
        STATICResponse = FALSE
      )
Statss <- .wilTest(
    ReferencesJ = SimulatedReferenceSpace[, 4],
    iplus = sample(
        1:nrow(SimulatedReferenceSpace),
        size = nrow(SimulatedReferenceSpace) / 20
    ),
    iminus = sample(
        1:nrow(SimulatedReferenceSpace),
        size = nrow(SimulatedReferenceSpace) / 20
    ),
    STATICResponse = TRUE
)

test_that("'.wilTest' has to return a single numeric value",{
  expect_is(loggedPs,"numeric")
  expect_length(loggedPs, 1)
  expect_is(Statss,"numeric")
  expect_length(Statss, 1)
})

test_that("'.wilTest' has to return minues-logged-pvalue",{
    expect_lte(2^-abs(loggedPs), 1)
    expect_gte(2^-abs(loggedPs), 0)
})

test_that("'.wilTest' returns normalized stats",{
    expect_lte(Statss, 1)
    expect_gte(Statss, -1)
})
