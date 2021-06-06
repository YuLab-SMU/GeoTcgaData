#' Convert  ENSEMBL gene id to gene Symbol in TCGA
#'
#' @param profiles a data.frame
#'
#' @return a data.frame, gene symbols and their expression value
#' @export
#'
#' @examples
#' result <- id_conversion(profile)
id_conversion<-function(profiles){  
    rownames(profiles) <- unlist(lapply(rownames(profiles), function(x) unlist(strsplit(x,"\\."))[1]))
    file3 <- hgnc
    rownames(file3)<-file3[,5]
    genes <- intersect(rownames(profiles), file3[,5])
    mat <- file3[genes, c(1,5)]
    profiles_new <- profiles[mat[,2],]
    rownames(profiles_new) <- mat[,1]
    profiles_new
}

