
## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1") utils::globalVariables("hgnc_file")

#' Gene id conversion types
#'
#' @return a vector 
#' @export
#'
#' @examples
#' id_ava()
id_ava <- function() {
    return(colnames(hgnc_file))
}


# id_conversion_one <- function(from,to,ID) {
    # if(from == "RefSeq_ID") from <- "refseq_accession"
    # if(from == "Ensembl_ID") from <- "ensembl_gene_id"
    # if(from == "NCBI_Gene_ID") from <- "NCBI Gene ID"
    # if(from == "UCSC_ID") from <- "ucsc_id"
    # if(from == "UniProt_ID") from <- "uniprot_ids"
    
    # if(to == "RefSeq_ID") to <- "refseq_accession"
    # if(to == "Ensembl_ID") to <- "ensembl_gene_id"
    # if(to == "NCBI_Gene_ID") to <- "NCBI Gene ID"
    # if(to == "UCSC_ID") to <- "ucsc_id"
    # if(to == "UniProt_ID") to <- "uniprot_ids"    

    # if(!(from %in% id_ava())) {
        # stop('from must be one of "id_ava()"')
    # }
    # if(!(to %in% id_ava())) {
        # stop('to must be one of "id_ava()"')
    # }
  
    
    # hgnc_filee <- hgnc_file
    # rownames(hgnc_filee) = hgnc_filee[,from]
    # if(ID != "" & ID %in% hgnc_filee[,from]) {
        # con <- hgnc_filee[ID,to]
    # } else { con <- "not available"}
    # return(con)
# }

#' Gene id conversion
#'
#' @param from one of "id_ava()"
#' @param to one of "id_ava()"
#' @param IDs the gene id which needed to convert
#'
#' @return a vector of genes
#' @export
#'
#' @examples
#' id_conversion_vector("symbol","Ensembl_ID",c("A2ML1","A2ML1-AS1","A4GALT","A12M1","AAAS"))
id_conversion_vector <- function(from,to,IDs) {
    if(from == "RefSeq_ID") from <- "refseq_accession"
    if(from == "Ensembl_ID") from <- "ensembl_gene_id"
    if(from == "NCBI_Gene_ID") from <- "NCBI Gene ID"
    if(from == "UCSC_ID") from <- "ucsc_id"
    if(from == "UniProt_ID") from <- "uniprot_ids"
    
    if(to == "RefSeq_ID") to <- "refseq_accession"
    if(to == "Ensembl_ID") to <- "ensembl_gene_id"
    if(to == "NCBI_Gene_ID") to <- "NCBI Gene ID"
    if(to == "UCSC_ID") to <- "ucsc_id"
    if(to == "UniProt_ID") to <- "uniprot_ids"    

    if(!(from %in% id_ava())) {
        stop('from must be one of "id_ava()"')
    }
    if(!(to %in% id_ava())) {
        stop('to must be one of "id_ava()"')
    }
    
    hgnc_filee <- hgnc_file
    rownames(hgnc_filee) = hgnc_filee[,from]
    results_id <- rep("not available", length(IDs))
    results_mat <- cbind(IDs, results_id)
    rownames(results_mat) <- IDs
    ID_hg <- intersect(IDs, hgnc_filee[,from])
    results_mat[ID_hg, 2] <- hgnc_filee[ID_hg, to]
    return(results_mat[, 2]) 
}

# id_conversion_vector <- function(from,to,IDs) {
  # results_id <-  IDs
  # for(i in seq_len(length(IDs))) {
    # results_id[i] <- id_conversion_one(from,to,IDs[i])
  # }
  # return(results_id)
# }
