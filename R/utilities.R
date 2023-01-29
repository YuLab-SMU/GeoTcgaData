#' Get Microarray matrix data from GEO
#'
#' @param gse GSE number, such as GSE781.
#'
#' @return a list of matrix
#' @export
#'
#' @examples
#' \donttest{
#' arraylist <- get_geo_array("GSE781")
#' }
get_geo_array <- function(gse) {
    gse <- GEOquery::getGEO(gse, GSEMatrix = FALSE, AnnotGPL = TRUE)
    gselist <- vector("list", length(GEOquery::GPLList(gse)))
    names(gselist) <- names(GEOquery::GPLList(gse))
    gsmplatforms <- lapply(GEOquery::GSMList(gse),
        function(x) {GEOquery::Meta(x)$platform_id})
    for (i in seq_len(length(gselist))) {
        gsmlist <- BiocGenerics::Filter(function(gsm) {
            GEOquery::Meta(gsm)$platform_id==names(gselist)[i]},
        GEOquery::GSMList(gse))
        probesets <- GEOquery::Table(GEOquery::GPLList(gse)[[1]])$ID
        data.matrix <- do.call('cbind',
            lapply(gsmlist, function(x) {tab <- GEOquery::Table(x)
                                        mymatch <- match(probesets,tab$ID_REF)
                                        return(tab$VALUE[mymatch])
                }
            )
        )
        data.matrix <- apply(data.matrix,2,
            function(x) {as.numeric(as.character(x))})
        gpl <- gse@gpls[names(gselist)[i]]
        gpl <- gpl[[1]]@dataTable@table
        genes <- gpl[match(probesets, gpl[, "ID"]), "Gene Symbol"]
        rownames(data.matrix) <- genes
        colnames(data.matrix) <- names(gsmlist)
        gselist[[i]] <- data.matrix
    }
}

#' Preprocess of Microarray data
#'
#' @param x matrix of Microarray data
#' @param missing_value Method to  impute missing expression data,
#' one of "zero" and "knn".
#' @param string a string, sep of the gene
#'
#' @return matrix
#' @export 
#'
#' @examples
#' \donttest{
#' arraylist <- get_geo_array("GSE781")
#' arraylist <- lapply(arraylist, array_preprocess)
#' }
array_preprocess <- function(x, missing_value = "knn", string = " /// ") {
    ## filter
    x <- x[!is.na(rownames(x)), ]
    x <- x[rownames(x) != "", ]
    aa <- rowSums(is.na(x))
    x <- x[aa < ncol(x)/2, ]

    ## impute missing
    if (missing_value == "zero") {
        x[is.na(x)] <- 0
    } else {
        x <- quiet(impute::impute.knn(x)$data)
    }

    ## log
    qx <- as.numeric(stats::quantile(x, 
        c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
    LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
    if (LogC) {
        x[which(x <= 0)] <- 0.0001
        x <- log2(x)
    }

    ## gene id conversion
    x <- cbind(rownames(x), x)
    x <- rep1(x, string)
    # gene_ave(x)
}

#' cluster probes of Microarray data
#'
#' @param x matrix of Microarray data, the first is the name of the gene, 
#' and the others are the expression value.
#' @param clusterCutoff Pearson correlation threshold 
#' to cut off the hierarchical tree.
#' @importFrom stats as.dist
#' @importFrom stats cor
#' @importFrom stats cutree
#' @importFrom stats hclust
#' @return data.frame
#' @export 
#'
#' @examples
#' \donttest{
#' arraylist <- get_geo_array("GSE781")
#' arraylist <- lapply(arraylist, array_preprocess)
#' arraylist_cluster <- lapply(arraylist, cluster_array)
#' }
cluster_array <- function(x, clusterCutoff = 0.7) {
    genes <- x[, 1]
    uniqueGenes <- unique(genes)
    x <- x[, -1]
    matlist <- vector("list", length(uniqueGenes))
    for (i in seq_len(length(uniqueGenes))) {
        gene <- uniqueGenes[i]
        probes <- which(genes == gene)
        mat <-  x[probes, ]
        if (length(probes) == 1) {
            rownames(mat) <- gene
            matlist[[i]] <- mat
        } else {
            probeCorrelation <- cor(t(mat),method = 'pearson')
            ClusterResults <- hclust(as.dist(1-probeCorrelation), 
                method = "complete", members = NULL)
            #plot(ClusterResults)
            Clusters <- cutree(ClusterResults, h = clusterCutoff)
            clusterDf <- matrix(0, length(unique(Clusters)), ncol(mat)) |> 
                as.data.frame()
            for (j in seq_len(length(unique(Clusters)))) {
                tmpGeneProbes <-  which(Clusters == j)
                if (length(tmpGeneProbes) > 1) {
                    clusterDf[j, ] <- colMeans(mat[tmpGeneProbes,])
                } else {
                    clusterDf[j, ] <- mat[tmpGeneProbes,]
                }        
            }
            if (nrow(clusterDf) > 1) {
                rownames(clusterDf) <- paste(gene, seq(nrow(clusterDf)), sep = "_")
            } else {
                rownames(clusterDf) <- gene
            }
            
            colnames(clusterDf) <- colnames(mat)    
            matlist[[i]] <- clusterDf
        }
    }
    matlist <- do.call("rbind", matlist)
}