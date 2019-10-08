countToTpm <- function(counts, effLen)
{
    rate <- log(counts) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
}
countToFpkm <- function(counts, effLen)
{
    N <- sum(counts)
    exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

fpkmToTpm <- function(fpkm)
{
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

countToEffCounts <- function(counts, len, effLen)
{
    counts * (len / effLen)
}

# if we have fpkm, then we can easily get the rate of counts/sum(counts).
# we can't get the real count value.
fpkmToCount <- function(fpkm, effLen, N=1e9)
{
	#rate <- (fpkm * effLen)/10^9
	rate <- exp(log(fpkm) + log(effLen) - log(1e9))
	counts <- rate * N
}




#' Convert count to FPKM
#'
#' @param counts_matrix a matrix, colnames of counts_matrix are sample name,
#' rownames of counts_matrix are gene symbols
#'
#' @return  a matrix
#' @export
#'
#' @examples
#' lung_squ_count2 <- matrix(c(1,2,3,4,5,6,7,8,9),ncol=3)
#' rownames(lung_squ_count2) <- c("DISC1","TCOF1","SPPL3")
#' colnames(lung_squ_count2) <- c("sample1","sample2","sample3")
#' jieguo <- countToFpkm_matrix(lung_squ_count2)
countToFpkm_matrix <- function(counts_matrix) {
    genes_count <- intersect(rownames(counts_matrix),gene_loc_len[,1])
	counts_matrix_new <- counts_matrix[genes_count,]
	gene_loc_len_new <- gene_loc_len[genes_count,]
	genes_length <- as.numeric(gene_loc_len_new[,4])
	counts_matrix_new2 <- counts_matrix_new
	for(i in seq_len(dim(counts_matrix_new2)[2])) {
	    counts_matrix_new2[,i] <- countToFpkm(as.numeric(counts_matrix_new2[,i]),genes_length)
	}
    return(counts_matrix_new2)
}


#' Convert count to Tpm
#'
#' @param counts_matrix a matrix, colnames of counts_matrix are sample name,
#' rownames of counts_matrix are gene symbols
#'
#' @return a matrix
#' @export
#'
#' @examples
#' lung_squ_count2 <- matrix(c(1,2,3,4,5,6,7,8,9),ncol=3)
#' rownames(lung_squ_count2) <- c("DISC1","TCOF1","SPPL3")
#' colnames(lung_squ_count2) <- c("sample1","sample2","sample3")
#' jieguo <- countToTpm_matrix(lung_squ_count2)
countToTpm_matrix <- function(counts_matrix) {
    genes_count <- intersect(rownames(counts_matrix),gene_loc_len[,1])
	counts_matrix_new <- counts_matrix[genes_count,]
	gene_loc_len_new <- gene_loc_len[genes_count,]
	genes_length <- as.numeric(gene_loc_len_new[,4])
	counts_matrix_new2 <- counts_matrix_new
	for(i in seq_len(dim(counts_matrix_new2)[2])) {
	    counts_matrix_new2[,i] <- countToTpm(as.numeric(counts_matrix_new2[,i]),genes_length)
	}
    return(counts_matrix_new2)
}


#' Convert fpkm to Tpm
#'
#' @param fpkm_matrix a matrix, colnames of fpkm_matrix are sample name,
#' rownames of fpkm_matrix are gene symbols
#'
#' @return a matrix
#' @export
#'
#' @examples
#' lung_squ_count2 <- matrix(c(0.11,0.22,0.43,0.14,0.875,0.66,0.77,0.18,0.29),ncol=3)
#' rownames(lung_squ_count2) <- c("DISC1","TCOF1","SPPL3")
#' colnames(lung_squ_count2) <- c("sample1","sample2","sample3")
#' jieguo <- countToTpm_matrix(lung_squ_count2)
fpkmToTpm_matrix <- function(fpkm_matrix) {
    genes_count <- intersect(rownames(fpkm_matrix),gene_loc_len[,1])
    fpkm_matrix_new <- fpkm_matrix[genes_count,]
    gene_loc_len_new <- gene_loc_len[genes_count,]
    genes_length <- as.numeric(gene_loc_len_new[,4])
    fpkm_matrix_new2 <- fpkm_matrix_new
    for(i in seq_len(dim(fpkm_matrix_new2)[2])) {
        fpkm_matrix_new2[,i] <- fpkmToTpm(as.numeric(fpkm_matrix_new2[,i]))
    }
    return(fpkm_matrix_new2)
}





