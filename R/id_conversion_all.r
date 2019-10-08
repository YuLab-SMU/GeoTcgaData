
id_conversion_one <- function(from,to,ID) {
    if(!(from %in% c("symbol","RefSeq_ID","Ensembl_ID","NCBI_Gene_ID","UCSC_ID","UniProt_ID"))) {
	    stop('from must be one of "symbol","RefSeq_ID",
		"Ensembl_ID","NCBI_Gene_ID","UCSC_ID","UniProt_ID"')
	}
	if(!(from %in% c("symbol","RefSeq_ID","Ensembl_ID","NCBI_Gene_ID","UCSC_ID","UniProt_ID"))) {
	    stop('to must be one of "symbol","RefSeq_ID",
	    "Ensembl_ID","NCBI_Gene_ID","UCSC_ID","UniProt_ID"')
	}
	hgnc_file = hgnc_complete_matrix
	rownames(hgnc_file) = hgnc_file[,from]
	if(ID %in% hgnc_file[,from]) {
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



