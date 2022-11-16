test_that("can parse example fpkmToTpm_matrix", {
  # skip_on_cran()
  lung_squ_count2 <- matrix(c(0.11, 0.22, 0.43, 0.14, 
    0.875, 0.66, 0.77, 0.18, 0.29), ncol = 3)
  rownames(lung_squ_count2) <- c("DISC1", "TCOF1", "SPPL3")
  colnames(lung_squ_count2) <- c("sample1", "sample2", "sample3")
  result <- fpkmToTpm_matrix(lung_squ_count2)
  expect_equal(dim(lung_squ_count2), dim(result))
})
