#' arrayDiff
#'
#' @param df data.frame of the omic data
#' @param group a vector, group of samples.
#' @param method one of "limma", "ttest", "wilcox",
#' @export
arrayDiff <- function(df, group, method = "limma") {
    method <- match.arg(method, c("limma", "ttest", "wilcox"))
    if (method == "limma") {
        result <- Diff_limma(df, group)
    } else {
        groups <- unique(group)
        which1 <- which(group == groups[1])
        which2 <- which(group == groups[2])
        pvalue <- rep(0, nrow(df))
        if (method == "ttest") {
            for (i in seq_len(length(pvalue))) {
                pvalue[i] <- stats::t.test(df[i, which1], df[i, which2])$p.value
                }
        } else {
            for (i in seq_len(length(pvalue))) {
                pvalue[i] <- stats::wilcox.test(df[i, which1], df[i, which2])$p.value
            }
        }
        result <- data.frame(gene = rownames(df), pvalue = pvalue)
    }
    return(df)

}