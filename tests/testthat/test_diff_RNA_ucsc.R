test_that("can parse example diff_RNA_ucsc", {
    # skip_on_bioc()
    df <- matrix(rnbinom(400, mu = 4, size = 10), 25, 16)
    df <- as.data.frame(df)
    rownames(df) <- paste0("gene", 1:25)
    colnames(df) <- paste0("sample", 1:16)
    df <- log2(df + 1)
    group <- sample(c("group1", "group2"), 16, replace = TRUE)
    df <- cbind(rownames(df), df)
    result <- diff_RNA_ucsc(ucsc = df, group = group, 
        filte = FALSE, method = "limma")
    expect_true( "P.Value"    %in% colnames(result))
})
