test_that("can parse example id_ava", {
  skip_on_cran()
  expect_equal(id_ava(), colnames(hgnc_file))
})
