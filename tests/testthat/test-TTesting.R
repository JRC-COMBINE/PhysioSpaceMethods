context("Test TTest wrapper function used by calculatePhysioMap")

SimulatedReferenceSpace <- matrix(rnorm(n = 10000, mean = 0, sd = 100),
                                  ncol = 10, dimnames = list(1:1000,11:20))

loggedPs <- .tTest(
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
Statss <- .tTest(
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

test_that("'.tTest' has to return a single numeric value",{
  expect_is(loggedPs,"numeric")
  expect_length(loggedPs, 1)
  expect_is(Statss,"numeric")
  expect_length(Statss, 1)
})

test_that("'.tTest' has to return minues-logged-pvalue",{
    expect_lte(2^-abs(loggedPs), 1)
    expect_gte(2^-abs(loggedPs), 0)
})

test_that("'.tTest' returns normalized stats",{
    expect_lte(Statss, 1)
    expect_gte(Statss, -1)
})

