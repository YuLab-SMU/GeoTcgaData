test_that("can parse example cal_mean_module", {
    data(geneExpress)
    data(module)
    result <- cal_mean_module(geneExpress, module)
    expect_equal( ncol(result) , ncol(geneExpress))
})