#' Convert  ENSEMBL gene id to gene Symbol in TCGA
#'
#' @param profiles a data.frame
#' @param toType one of 'keytypes(org.Hs.eg.db)'
#'
#' @return a data.frame, gene symbols and their expression value
#' @export
#'
#' @examples
#' library(org.Hs.eg.db)
#' profile <- GeoTcgaData::profile
#' result <- id_conversion(profile)
id_conversion<-function(profiles, toType = "SYMBOL"){  
    # rownames(profiles) <- unlist(lapply(rownames(profiles), function(x) unlist(strsplit(x,"\\."))[1]))
    # file3 <- hgnc
    # rownames(file3)<-file3[,5]
    # genes <- intersect(rownames(profiles), file3[,5])
    # mat <- file3[genes, c(1,5)]
    # profiles_new <- profiles[mat[,2],]
    # rownames(profiles_new) <- mat[,1]
    # profiles_new
    ## use clusterProfiler::bitr to convert gene id
    rownames(profiles) <- unlist(lapply(rownames(profiles), function(x) unlist(strsplit(x,"\\."))[1]))
    genes <- clusterProfiler::bitr(rownames(profiles), fromType = "ENSEMBL", 
        toType = toType, OrgDb = org.Hs.eg.db::org.Hs.eg.db, drop = FALSE)[, 2]
    profiles_new <- profiles[!is.na(genes), ]
    rownames(profiles_new) <- genes[!is.na(genes)]
    profiles_new
}

