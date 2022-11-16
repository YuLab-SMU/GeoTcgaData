test_that("can parse example differential_cnv", {
  # skip_on_cran()
  cnv <- matrix(c(
    -1.09150, -1.47120, -0.87050, -0.50880,
    -0.50880, 2.0, 2.0, 2.0, 2.0, 2.0, 2.601962, 2.621332, 2.621332,
    2.621332, 2.621332, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
    2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0
  ), nrow = 5)
  cnv <- as.data.frame(cnv)
  rownames(cnv) <- c("AJAP1", "FHAD1", "CLCNKB", "CROCCP2", "AL137798.3")
  colnames(cnv) <- c(
    "TCGA-DD-A4NS-10A-01D-A30U-01", "TCGA-ED-A82E-01A-11D-A34Y-01",
    "TCGA-WQ-A9G7-01A-11D-A36W-01", "TCGA-DD-AADN-01A-11D-A40Q-01",
    "TCGA-ZS-A9CD-10A-01D-A36Z-01", "TCGA-DD-A1EB-11A-11D-A12Y-01"
  )
  cnv_chi_file <- prepare_chi(cnv)
  chiResult <- differential_cnv(cnv_chi_file)
  expect_true("P.Value" %in% colnames(chiResult))
})
