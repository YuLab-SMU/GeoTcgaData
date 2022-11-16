test_that("can parse example Merge_methy_tcga", {
  # skip_on_cran()
  merge_result <- Merge_methy_tcga(system.file(file.path("extdata", "methy"), 
    package = "GeoTcgaData"))
  expect_equal(names(merge_result), c("methyResult", "cpg_info"))
})

test_that("can parse example Diff_limma", {
  # skip_on_cran()
  df <- matrix(runif(200), 25, 8)
  df <- as.data.frame(df)
  rownames(df) <- paste0("gene", 1:25)
  colnames(df) <- paste0("sample", 1:8)
  group <- sample(c("group1", "group2"), 8, replace = TRUE)
  result <- Diff_limma(df = df, group = group)
  expect_true( "P.Value"  %in% colnames(result))
})
