test_that("can parse example diff_CNV", {
    # skip_on_bioc()
    aa <- matrix(sample(c(0, 1, -1), 200, replace = TRUE), 25, 8)
    rownames(aa) <- paste0("gene", 1:25)
    colnames(aa) <- paste0("sample", 1:8)
    sampleGroup <- sample(c("A", "B"), ncol(aa), replace = TRUE)
    result <- diff_CNV(aa, sampleGroup)
    expect_true( "P.Value"    %in% colnames(result))
})
