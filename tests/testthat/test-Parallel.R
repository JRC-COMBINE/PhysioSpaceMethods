context("Test initializing parallel processing")

NumOfCoresToTest <- parallel::detectCores()-1
CL <- parallelInitializer(NumbrOfCores = NumOfCoresToTest)

test_that("'parallelInitializer' output test",{
  expect_is(CL,"SOCKcluster")
  expect_length(CL, NumOfCoresToTest)
  expect_error(parallelInitializer(NumbrOfCores =
                                       parallel::detectCores()+1))
})
