context("Test plotting a heatmap")

randMatInpt <-
    matrix(data = rnorm(n = 4000, mean = 10, sd = 20), nrow = 400)
rownames(randMatInpt) <- paste("ROWS", 1:400)
colnames(randMatInpt) <- paste("Sample", 1:10)

randMatRef <-
    matrix(data = rnorm(n = 12000, mean = 10, sd = 20), nrow = 400)
rownames(randMatRef) <- paste("ROWS", 1:400)
colnames(randMatRef) <- paste("Space", 1:30)

res <-
    calculatePhysioMap(InputData = randMatInpt, Space = randMatRef)



test_that("'PhysioHeatmap' testing...",{
    expect_true(PhysioHeatmap(PhysioResults = res,
                                 main = "Heatmap Testing"))
    expect_true(PhysioHeatmap(
        PhysioResults = res,
        main = "Heatmap Testing",
        ColorLevels = 3
    ))
    expect_true(PhysioHeatmap(
        PhysioResults = res,
        main = "Heatmap Testing",
        SpaceClustering = TRUE,
        Space = randMatRef
    ))
    expect_true(PhysioHeatmap(
        PhysioResults = res,
        main = "Heatmap Testing",
        ReducedPlotting = 2
    ))
})
