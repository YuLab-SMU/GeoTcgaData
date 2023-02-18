test_that("can parse example differential_array", {
    # skip_on_bioc()
    arrayData <- matrix(runif(200), 25, 8)
    rownames(arrayData) <- paste0("gene", 1:25)
    colnames(arrayData) <- paste0("sample", 1:8)
    group <- c(rep("group1", 4), rep("group2", 4))
    result <- differential_array(df = arrayData, group = group)
    expect_true( "P.Value"    %in% colnames(result))
})
