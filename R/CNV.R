#' Do difference analysis of gene level copy number variation data
#'
#' @param cnvData data.frame of CNV data
#' @param sampleGroup vector of sample group
#' @param ... parameters for fisher.test
#' @export
#'
#' @examples
#' \dontrun{
# use TCGAbiolinks to download TCGA data    
#' library(TCGAbiolinks)
#' query <- GDCquery(project = "TCGA-LGG",
#'                   data.category = "Copy Number Variation",
#'                   data.type = "Gene Level Copy Number Scores")
#' 
#' GDCdownload(query, method = "api", files.per.chunk = 5, directory = Your_Path)
#' data <- GDCprepare(query = query, 
#'                    save = TRUE, 
#'                    directory =  "Your_Path") 
#' 
#' class(data) <- "data.frame"
#' cnvData <- data[, -c(1,2,3)]
#' rownames(cnvData) <- data[, 1]
#' sampleGroup  = sample(c("A","B"), ncol(cnvData), replace = TRUE)
#' diffCnv <- diff_CNV(cnvData, sampleGroup)
#' }
diff_CNV <- function(cnvData, sampleGroup, ...) {
    type1 <- which(sampleGroup == names(table(sampleGroup))[1])
    type2 <- which(sampleGroup == names(table(sampleGroup))[2])
    pvalue <- rep(0, nrow(cnvData))
    estimate <- rep(0, nrow(cnvData))
    for (i in seq_len(nrow(cnvData))) {
        type1_freq <- table(as.character(cnvData[i, type1]))
        type2_freq <- table(as.character(cnvData[i, type2]))
        df <- data.frame(type1 = as.numeric(type1_freq[c("-1", "0", "1")]),
                         type2 = as.numeric(type2_freq[c("-1", "0", "1")]))
        df[is.na(df)] <- 0
        # rownames(df) <- c("-1", "0", "1")
        # df[2, ] <- df[3, ] + df[2, ]
        # df <- df[-3, ]
        fish <- stats::fisher.test(df, ...)
        pvalue[i] <- fish$p.value
        estimate[i] <- fish$estimate
    }
    names(pvalue) <- names(estimate) <- gsub("\\..*", "", rownames(cnvData))
    return(data.frame(pvalue, estimate))
}

diff_CNV_segment <- function(cnvData, sampleType) {
    
}





