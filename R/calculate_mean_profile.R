
#' Average the values of same genes in gene expression profile
#'
#' @param file_gene_ave a data.frame
#' @param k a number
#'
#' @return a data.frame, the values of same genes in gene expression profile
#' @export
#'
#' @examples
#' aa <- c("MARCH1", "MARC1", "MARCH1", "MARCH1", "MARCH1")
#' bb <- c(2.969058399, 4.722410064, 8.165514853, 8.24243893, 8.60815086)
#' cc <- c(3.969058399, 5.722410064, 7.165514853, 6.24243893, 7.60815086)
#' file_gene_ave <- data.frame(aa = aa, bb = bb, cc = cc)
#' colnames(file_gene_ave) <- c("Gene", "GSM1629982", "GSM1629983")
#'
#' result <- gene_ave(file_gene_ave, 1)
gene_ave <- function(file_gene_ave, k = 1) {
    x <- file_gene_ave[, -k]
    file_gene_ave <- as.matrix(file_gene_ave)
    rownames(file_gene_ave) <- file_gene_ave[, k]
    # x <- file_gene_ave
    ID <- rownames(file_gene_ave)
    ID <- factor(ID, levels = unique(ID))

    y <- rowsum(x, ID, reorder = FALSE, na.rm = TRUE)
    n <- rowsum(1L - is.na(x), ID, reorder = FALSE)
    return(y / n)
}
