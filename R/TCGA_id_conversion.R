#' Convert ENSEMBL gene id to gene Symbol in TCGA
#'
#' @param profiles a data.frame of gene expression data, 
#' each column is a sample, 
#' and each row is a gene. 
#' @param toType one of 'keytypes(org.Hs.eg.db)'
#'
#' @return a data.frame, gene symbols and their expression value
#' @export
#'
#' @examples
#' library(org.Hs.eg.db)
#' data(profile)
#' result <- id_conversion_TCGA(profile)
id_conversion_TCGA <- function(profiles, toType = "SYMBOL") {
    rownames(profiles) <- gsub("\\..*", "", rownames(profiles))
    genes <- clusterProfiler::bitr(rownames(profiles),
        fromType = "ENSEMBL",
        toType = toType, OrgDb = org.Hs.eg.db::org.Hs.eg.db, drop = FALSE
    )

    genes <- genes[!duplicated(genes[, 1]), ]
    rownames(genes) <- genes[, 1]
    profiles2 <- as.matrix(profiles)
    rownames(profiles2) <- genes[rownames(profiles), 2]
    return(profiles2)
}
