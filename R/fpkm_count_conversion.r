countToTpm_internal <- function(counts, effLen) {
    rate <- log(counts) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
}
countToFpkm_internal <- function(counts, effLen) {
    N <- sum(counts)
    exp(log(counts) + log(1e9) - log(effLen) - log(N))
}

fpkmToTpm_internal <- function(fpkm) {
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

countToEffCounts_internal <- function(counts, len, effLen) {
    counts * (len / effLen)
}

# if we have fpkm, then we can easily get the rate of counts/sum(counts).
# we can't get the real count value.
fpkmToCount_internal <- function(fpkm, effLen, N = 1e9) {
    # rate <- (fpkm * effLen)/10^9
    rate <- exp(log(fpkm) + log(effLen) - log(1e9))
    counts <- rate * N
}




#' Convert count to FPKM
#'
#' @param counts_matrix a matrix, colnames of counts_matrix are sample name,
#' rownames of counts_matrix are gene symbols
#' @param keyType keyType, one of keytypes(org.Hs.eg.db).
#' @param gene_cov data.frame of two column, the first column is gene length, 
#' the second column is gene GC content
#'
#' @return a matrix
#' @export
#'
#' @examples
#' data(gene_cov)
#' lung_squ_count2 <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), ncol = 3)
#' rownames(lung_squ_count2) <- c("DISC1", "TCOF1", "SPPL3")
#' colnames(lung_squ_count2) <- c("sample1", "sample2", "sample3")
#' result <- countToFpkm(lung_squ_count2,
#'     keyType = "SYMBOL",
#'     gene_cov = gene_cov
#' )
countToFpkm <- function(counts_matrix, keyType = "SYMBOL", gene_cov) {
    gene_cov2 <- gene_cov
    if (keyType != "ENTREZID") {
        genes_bitr <- clusterProfiler::bitr(rownames(gene_cov),
            fromType = "ENTREZID", toType = keyType,
            OrgDb = org.Hs.eg.db::org.Hs.eg.db, drop = TRUE
        )
        genes_bitr <- genes_bitr[!duplicated(genes_bitr[, 2]), ]
        gene_cov2 <- gene_cov[genes_bitr$ENTREZID, ]
        rownames(gene_cov2) <- genes_bitr[, 2]
    }
    genes_count <- intersect(rownames(counts_matrix), rownames(gene_cov2))
    counts_matrix_new <- counts_matrix[genes_count, ]
    gene_loc_len_new <- gene_cov2[genes_count, ]
    genes_length <- as.numeric(gene_loc_len_new[, 1])
    counts_matrix_new2 <- counts_matrix_new
    for (i in seq_len(dim(counts_matrix_new2)[2])) {
        counts_matrix_new2[, i] <- countToFpkm_internal(
            as.numeric(counts_matrix_new2[, i]),
            genes_length)
    }
    return(counts_matrix_new2)
}


#' Convert count to Tpm
#'
#' @param counts_matrix a matrix, colnames of counts_matrix are sample name,
#' rownames of counts_matrix are gene symbols
#' @param keyType keyType, one of keytypes(org.Hs.eg.db).
#' @param gene_cov data.frame of two column, the first column is gene length, 
#' the second column is gene GC content
#'
#' @return a matrix
#' @export
#'
#' @examples
#' data(gene_cov)
#' lung_squ_count2 <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), ncol = 3)
#' rownames(lung_squ_count2) <- c("DISC1", "TCOF1", "SPPL3")
#' colnames(lung_squ_count2) <- c("sample1", "sample2", "sample3")
#' result <- countToTpm(lung_squ_count2,
#'     keyType = "SYMBOL",
#'     gene_cov = gene_cov
#' )
countToTpm <- function(counts_matrix, keyType = "SYMBOL", gene_cov) {
    gene_cov2 <- gene_cov
    if (keyType != "ENTREZID") {
        genes_bitr <- clusterProfiler::bitr(rownames(gene_cov),
            fromType = "ENTREZID", toType = keyType,
            OrgDb = org.Hs.eg.db::org.Hs.eg.db, drop = TRUE
        )
        genes_bitr <- genes_bitr[!duplicated(genes_bitr[, 2]), ]
        gene_cov2 <- gene_cov[genes_bitr$ENTREZID, ]
        rownames(gene_cov2) <- genes_bitr[, 2]
    }
    genes_count <- intersect(rownames(counts_matrix), rownames(gene_cov2))
    counts_matrix_new <- counts_matrix[genes_count, ]
    gene_loc_len_new <- gene_cov2[genes_count, ]
    genes_length <- as.numeric(gene_loc_len_new[, 1])
    counts_matrix_new2 <- counts_matrix_new
    for (i in seq_len(dim(counts_matrix_new2)[2])) {
        counts_matrix_new2[, i] <- countToTpm_internal(
            as.numeric(counts_matrix_new2[, i]),
            genes_length)
    }
    return(counts_matrix_new2)
}


#' Convert fpkm to Tpm
#'
#' @param fpkm_matrix a matrix, colnames of fpkm_matrix are sample name,
#' rownames of fpkm_matrix are genes
#'
#' @return a matrix
#' @export
#'
#' @examples
#' lung_squ_count2 <- matrix(c(0.11, 0.22, 0.43, 0.14, 0.875,
#'     0.66, 0.77, 0.18, 0.29), ncol = 3)
#' rownames(lung_squ_count2) <- c("DISC1", "TCOF1", "SPPL3")
#' colnames(lung_squ_count2) <- c("sample1", "sample2", "sample3")
#' result <- fpkmToTpm(lung_squ_count2)
fpkmToTpm <- function(fpkm_matrix) {
    fpkm_matrix_new <- apply(fpkm_matrix, 2, fpkmToTpm_internal)
}
