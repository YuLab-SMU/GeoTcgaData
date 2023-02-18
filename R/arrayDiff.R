#' differential_array
#'
#' @param df data.frame of the omic data
#' @param group a vector, group of samples.
#' @param method one of "limma", "ttest", "wilcox"
#' @param adjust.method adjust.method.
#' @return data.frame
#' @export
#'
#' @examples
#' \donttest{
#' library(GeoTcgaData)
#' library(data.table)
#' # Use real GEO data as example
#' arrayData <- read.table("GSE54807_series_matrix.txt.gz",
#'     sep = "\t", header = TRUE,
#'         fill=TRUE, comment.char = "!", check.names=FALSE)
#' gpl <- fread("GPL6244-17930.txt", sep = "\t", header = TRUE)
#' gpl <- gpl[, c("ID", "gene_assignment")]
#' class(gpl) <- "data.frame"
#'
#' for (i in seq_len(nrow(gpl))) {
#'         aa <- strsplit(gpl[i, 2], " // ")[[1]][5]
#'         gpl[i, 2] <- as.character(strsplit(aa, " /// ")[[1]][1])
#' }
#' gpl[,1] <- as.character(gpl[,1])
#' arrayData[, 1] <- as.character(arrayData[, 1])
#' rownames(gpl) <- gpl[, 1]
#' arrayData[, 1] <- gpl[arrayData[, 1], 2]
#'
#'
#' arrayData <- repRemove(arrayData," /// ")
#'
#' # Remove rows that do not correspond to genes
#' arrayData <- arrayData[!is.na(arrayData[, 1]), ]
#' arrayData <- arrayData[!arrayData[, 1] == "", ]
#' arrayData <- arrayData[!arrayData[, 1] == "---", ]
#'
#'
#' arrayData <- arrayData[order(arrayData[, 1]), ]
#' arrayData <- gene_ave(arrayData, 1)
#'
#' keep <- apply(arrayData, 1, function(x) sum(x < 1) < (length(x)/2))
#' arrayData <- arrayData[keep, ]
#'
#' group <- c(rep("group1", 12), rep("group2", 12))
#' result <- differential_array(df = arrayData, group = group)
#' }
#' # Use random data as example
#' arrayData <- matrix(runif(200), 25, 8)
#' rownames(arrayData) <- paste0("gene", 1:25)
#' colnames(arrayData) <- paste0("sample", 1:8)
#' group <- c(rep("group1", 4), rep("group2", 4))
#' result <- differential_array(df = arrayData, group = group)
differential_array <- function(df, group, method = "limma", adjust.method = "BH") {
    method <- match.arg(method, c("limma", "ttest", "wilcox"))
    if (method == "limma") {
        result <- differential_limma(df, group, adjust.method = adjust.method)
    } else {
        groups <- unique(group)
        which1 <- which(group == groups[1])
        which2 <- which(group == groups[2])
        P.Value <- rep(0, nrow(df))
        if (method == "ttest") {
            for (i in seq_len(length(P.Value))) {
                P.Value[i] <- stats::t.test(df[i, which1],
                    df[i, which2])$p.value
            }
        } else {
            for (i in seq_len(length(P.Value))) {
                P.Value[i] <- stats::wilcox.test(as.numeric(df[i, which1]),
                    as.numeric(df[i, which2]))$p.value
            }
        }
        adj.P.Val <- stats::p.adjust(P.Value, method = adjust.method)
        result <- data.frame(gene = rownames(df),
            P.Value = P.Value, adj.P.Val = adj.P.Val)
    }
    return(result)
}
