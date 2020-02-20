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


id_conversion_one <- function(from,to,ID) {
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
  
    
   
    rownames(hgnc_file) = hgnc_file[,from]
    if(ID != "" & ID %in% hgnc_file[,from]) {
        con <- hgnc_file[ID,to]
    } else { con <- "not available"}
    return(con)
}

#' Gene id conversion
#'
#' @param from one of "symbol","RefSeq_ID","Ensembl_ID","NCBI_Gene_ID","UCSC_ID","UniProt_ID"
#' @param to one of "symbol","RefSeq_ID","Ensembl_ID","NCBI_Gene_ID","UCSC_ID","UniProt_ID"
#' @param IDs the gene id which needed to convert
#'
#' @return a vector of genes
#' @export
#'
#' @examples
#' id_conversion_vector("symbol","Ensembl_ID",c("A2ML1","A2ML1-AS1","A4GALT","A12M1","AAAS"))
id_conversion_vector <- function(from,to,IDs) {
  results_id <-  IDs
  for(i in seq_len(length(IDs))) {
    results_id[i] <- id_conversion_one(from,to,IDs[i])
  }
  return(results_id)
}



