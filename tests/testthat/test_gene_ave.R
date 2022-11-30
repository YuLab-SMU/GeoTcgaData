test_that("can parse example gene_ave", {
    # skip_on_bioc()
    aa <- c("MARCH1", "MARC1", "MARCH1", "MARCH1", "MARCH1")
    bb <- c(2.969058399, 4.722410064, 8.165514853, 8.24243893, 8.60815086)
    cc <- c(3.969058399, 5.722410064, 7.165514853, 6.24243893, 7.60815086)
    file_gene_ave <- data.frame(aa = aa, bb = bb, cc = cc)
    colnames(file_gene_ave) <- c("Gene", "GSM1629982", "GSM1629983")
    result <- gene_ave(file_gene_ave, 1)
    expect_equal(sort(unique(file_gene_ave[, 1])),    
        sort(unique(rownames(result))))
})
