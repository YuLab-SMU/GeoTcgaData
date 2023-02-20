test_that("can parse example Merge_methy_tcga", {
    merge_result <- Merge_methy_tcga(system.file(file.path("extdata", "methy"), 
      package = "GeoTcgaData"))
    expect_equal(names(merge_result), c("methyResult", "cpg_info"))
})