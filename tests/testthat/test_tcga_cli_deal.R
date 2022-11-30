test_that("can parse example tcga_cli_deal", {
    # skip_on_bioc()
    result <- tcga_cli_deal(system.file(file.path("extdata", "tcga_cli"), 
        package = "GeoTcgaData"))
    expect_true("Dead"    %in% result[, 3])
})
