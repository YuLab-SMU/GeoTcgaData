#' Do difference analysis of gene level copy number variation data
#'
#' @param cnvData data.frame of CNV data
#' @param sampleGroup vector of sample group
#' @param method one of "Chisquare", "fisher",
#' and "CATT"(Cochran-Armitage trend test)
#' @param adjust.method adjust.method.
#' @param ... parameters for "Chisquare", "fisher",
#' and "CATT"(Cochran-Armitage trend test)
#' @return data.frame with pvalue and estimate
#' @export
#'
#' @examples
#' \donttest{
#' # use TCGAbiolinks data as example
#' library(TCGAbiolinks)
#' query <- GDCquery(
#'         project = "TCGA-ACC",
#'         data.category = "Copy Number Variation",
#'         data.type = "Gene Level Copy Number",
#'         access = "open"
#' )
#' GDCdownload(query)
#' cnvData <- GDCprepare(query)
#' aa <- assays(cnvData)$copy_number
#' bb <- aa
#' aa[bb == 2] <- 0
#' aa[bb < 2] <- -1
#' aa[bb > 2] <- 1
#' sampleGroup <- sample(c("A", "B"), ncol(cnvData), replace = TRUE)
#' diffCnv <- diff_CNV(aa, sampleGroup)
#'
#' # Use sangerbox CNV data as example
#' cnvData <- fread("Merge_GeneLevelCopyNumber.txt")
#' class(cnvData) <- "data.frame"
#' rownames(cnvData) <- cnvData[, 1]
#' cnvData <- cnvData[, -c(1, 2, 3)]
#' sampleGroup <- sample(c("A", "B"), ncol(cnvData), replace = TRUE)
#' diffCnv <- diff_CNV(cnvData, sampleGroup)
#' }
#' # use random data as example
#' aa <- matrix(sample(c(0, 1, -1), 200, replace = TRUE), 25, 8)
#' rownames(aa) <- paste0("gene", 1:25)
#' colnames(aa) <- paste0("sample", 1:8)
#' sampleGroup <- sample(c("A", "B"), ncol(aa), replace = TRUE)
#' diffCnv <- diff_CNV(aa, sampleGroup)

diff_CNV <- function(cnvData, sampleGroup,
                    method = "Chisquare",
                    adjust.method = "BH", ...) {
    type1 <- which(sampleGroup == unique(sampleGroup)[1])
    type2 <- which(sampleGroup == unique(sampleGroup)[2])
    pvalue <- rep(1, nrow(cnvData))
    estimate <- rep(0, nrow(cnvData))
    for (i in seq_len(nrow(cnvData))) {
        type1_freq <- table(as.character(cnvData[i, type1]))
        type2_freq <- table(as.character(cnvData[i, type2]))
        df <- data.frame(
            type1 = as.numeric(type1_freq[c("-1", "0", "1")]),
            type2 = as.numeric(type2_freq[c("-1", "0", "1")])
        )
        df[is.na(df)] <- 0
        df <- df[rowSums(df) > 0, ]
        if (nrow(df) > 2) {
            if (method == "fisher") {
                fish <- stats::fisher.test(df, ...)
                pvalue[i] <- fish$p.value
                if (nrow(df) == 2) {
                        estimate[i] <- fish$estimate
                }
            }

            if (method == "Chisquare") {
                    pvalue[i] <- stats::chisq.test(df, ...)$p.value
            }

            if(method == "CATT") {
                    pvalue[i] <- CATT::CATT(table = t(df), ...)$p.value
            }
        }
    }

    #     fish <- tryCatch({
    #         stats::fisher.test(df, ...)
    #     }, warning = function(w){
    #     }, error = function(e){
    #             list(p.value = 1, estimate = 0)
    #     },finally = {
    #     })
    #     P.Value[i] <- fish$p.value
    # }
    adj.P.Val <- stats::p.adjust(pvalue, method = adjust.method)
    gene <- gsub("\\..*", "", rownames(cnvData))
    result <- data.frame(gene = gene, P.Value = pvalue,
        adj.P.Val = adj.P.Val, estimate = estimate)
    rownames(result) <- gene
    result
}

# diff_CNV_segment <- function(cnvData, sampleType) {

# }
