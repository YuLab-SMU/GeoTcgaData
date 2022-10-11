#' Preparer file for chi-square test
#'
#' @param cnv result of ann_merge()
#'
#' @return a matrix
#' @export
#'
#' @examples
#' cnv <- matrix(c(
#'   -1.09150, -1.47120, -0.87050, -0.50880,
#'   -0.50880, 2.0, 2.0, 2.0, 2.0, 2.0, 2.601962, 2.621332, 2.621332,
#'   2.621332, 2.621332, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
#'   2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0
#' ), nrow = 5)
#' cnv <- as.data.frame(cnv)
#' rownames(cnv) <- c("AJAP1", "FHAD1", "CLCNKB", "CROCCP2", "AL137798.3")
#' colnames(cnv) <- c(
#'   "TCGA-DD-A4NS-10A-01D-A30U-01", "TCGA-ED-A82E-01A-11D-A34Y-01",
#'   "TCGA-WQ-A9G7-01A-11D-A36W-01", "TCGA-DD-AADN-01A-11D-A40Q-01",
#'   "TCGA-ZS-A9CD-10A-01D-A36Z-01", "TCGA-DD-A1EB-11A-11D-A12Y-01"
#' )
#' cnv_chi_file <- prepare_chi(cnv)
prepare_chi <- function(cnv) {
  file1 <- cnv
  sampless <- rep(0, ncol(file1))
  for (i in seq_len(length(sampless)))
  {
    a <- unlist(strsplit(colnames(file1)[i], "-"))[4]
    sampless[i] <- substring(a, 1, 1)
  }

  cnv_chi <- matrix(0, nrow(cnv), 4)
  rownames(cnv_chi) <- rownames(cnv)
  colnames(cnv_chi) <- c("normalCNV", "normalWild", "tumorCNV", "tumorWild")
  for (i in seq_len(nrow(file1)))
  {
    normalCNV <- 0
    normalWild <- 0
    tumorCNV <- 0
    tumorWild <- 0
    for (j in seq_len(ncol(file1)))
    {
      if ((sampless[j] == "1") && (abs(as.numeric(file1[i, j]) - 2) > 0.5)) {
        normalCNV <- normalCNV + 1
      }
      if ((sampless[j] == "1") && (abs(as.numeric(file1[i, j]) - 2) <= 0.5)) {
        normalWild <- normalWild + 1
      }
      if ((sampless[j] == "0") && (abs(as.numeric(file1[i, j]) - 2) > 0.5)) {
        tumorCNV <- tumorCNV + 1
      }
      if ((sampless[j] == "0") && (abs(as.numeric(file1[i, j]) - 2) <= 0.5)) {
        tumorWild <- tumorWild + 1
      }
    }
    cnv_chi[i, ] <- c(normalCNV, normalWild, tumorCNV, tumorWild)
  }
  return(cnv_chi)
}


#' Do chi-square test to find differential genes
#'
#' @param rt result of prepare_chi()
#'
#' @return a matrix
#' @export
#'
#' @examples
#' result3 <- matrix(c(
#'   -1.09150, -1.47120, -0.87050, -0.50880,
#'   -0.50880, 2.0, 2.0, 2.0, 2.0, 2.0, 2.601962, 2.621332, 2.621332,
#'   2.621332, 2.621332, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
#'   2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0
#' ), nrow = 5)
#' rownames(result3) <- c("AJAP1", "FHAD1", "CLCNKB", "CROCCP2", "AL137798.3")
#' colnames(result3) <- c(
#'   "TCGA-DD-A4NS-10A-01D-A30U-01", "TCGA-ED-A82E-01A-11D-A34Y-01",
#'   "TCGA-WQ-A9G7-01A-11D-A36W-01", "TCGA-DD-AADN-01A-11D-A40Q-01",
#'   "TCGA-ZS-A9CD-10A-01D-A36Z-01", "TCGA-DD-A1EB-11A-11D-A12Y-01"
#' )
#' rt <- prepare_chi(result3)
#' chiResult <- differential_cnv(rt)
differential_cnv <- function(rt) {
  outTab <- matrix(0, nrow(rt), 4)
  colnames(outTab) <- c("normalRatio", "tumorRatio", "Stat", "P.Value")
  rownames(outTab) <- rownames(rt)
  for (i in seq_len(nrow(rt))) {
    x <- matrix(rt[i, ], ncol = 2)
    chiTest <- stats::chisq.test(x)
    normalRatio <- rt[i, 1] / (rt[i, 1] + rt[i, 2])
    tumorRatio <- rt[i, 3] / (rt[i, 3] + rt[i, 4])
    Gene <- row.names(rt[i, ])
    Stat <- chiTest$statistic
    P.Value <- chiTest$p.value
    outTab[i, ] <- c(normalRatio, tumorRatio, Stat, P.Value)
  }
  outTab <- as.data.frame(outTab)
  outTab$adj.P.Val <- stats::p.adjust(outTab$P.Value, method = "bonferroni")
  outTab
}
