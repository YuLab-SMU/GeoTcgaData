test_that("can parse example cal_mean_module", {
  skip_on_cran()
  result <- cal_mean_module(geneExpress, module)
  expect_equal( ncol(result) , ncol(geneExpress))
})