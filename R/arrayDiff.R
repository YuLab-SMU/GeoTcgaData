#' Differential analysis of Microarray data
#'
#' @param df data.frame of the omic data, each column is a sample, 
#' and each row is a gene. 
#' @param group a vector, group of samples.
#' @param method method to do differential analysis, 
#' one of "limma", "ttest", "wilcox".
#' @param adjust.method adjust.method, one of "holm", "hochberg", "hommel", 
#' "bonferroni", "BH", "BY", "fdr", and "none". 
#' @return data.frame
#' @export
#'
#' @examples
#' \donttest{
#' library(GeoTcgaData)
#' library(data.table)
#' # Use real GEO data as example
#' arrayData <- read.table("GSE54807_series_matrix.txt.gz",
#'     sep = "\t", header = TRUE,
#'         fill=TRUE, comment.char = "!", check.names=FALSE)
#' gpl <- fread("GPL6244-17930.txt", sep = "\t", header = TRUE)
#' gpl <- gpl[, c("ID", "gene_assignment")]
#' class(gpl) <- "data.frame"
#'
#' for (i in seq_len(nrow(gpl))) {
#'         aa <- strsplit(gpl[i, 2], " // ")[[1]][5]
#'         gpl[i, 2] <- as.character(strsplit(aa, " /// ")[[1]][1])
#' }
#' gpl[,1] <- as.character(gpl[,1])
#' arrayData[, 1] <- as.character(arrayData[, 1])
#' rownames(gpl) <- gpl[, 1]
#' arrayData[, 1] <- gpl[arrayData[, 1], 2]
#'
#'
#' arrayData <- repRemove(arrayData," /// ")
#'
#' # Remove rows that do not correspond to genes
#' arrayData <- arrayData[!is.na(arrayData[, 1]), ]
#' arrayData <- arrayData[!arrayData[, 1] == "", ]
#' arrayData <- arrayData[!arrayData[, 1] == "---", ]
#'
#'
#' arrayData <- arrayData[order(arrayData[, 1]), ]
#' arrayData <- gene_ave(arrayData, 1)
#'
#' keep <- apply(arrayData, 1, function(x) sum(x < 1) < (length(x)/2))
#' arrayData <- arrayData[keep, ]
#'
#' group <- c(rep("group1", 12), rep("group2", 12))
#' result <- differential_array(df = arrayData, group = group)
#' }
#' # Use random data as example
#' arrayData <- matrix(runif(200), 25, 8)
#' rownames(arrayData) <- paste0("gene", 1:25)
#' colnames(arrayData) <- paste0("sample", 1:8)
#' group <- c(rep("group1", 4), rep("group2", 4))
#' result <- differential_array(df = arrayData, group = group)
differential_array <- function(df, group, method = "limma", 
                                adjust.method = "BH") {
    method <- match.arg(method, c("limma", "ttest", "wilcox"))
    if (method == "limma") {
        result <- differential_limma(df, group, adjust.method = adjust.method)
    } else {
        groups <- unique(group)
        which1 <- which(group == groups[1])
        which2 <- which(group == groups[2])
        P.Value <- rep(0, nrow(df))
        if (method == "ttest") {
            for (i in seq_len(length(P.Value))) {
                P.Value[i] <- stats::t.test(df[i, which1],
                    df[i, which2])$p.value
            }
        } else {
            for (i in seq_len(length(P.Value))) {
                P.Value[i] <- stats::wilcox.test(as.numeric(df[i, which1]),
                    as.numeric(df[i, which2]))$p.value
            }
        }
        adj.P.Val <- stats::p.adjust(P.Value, method = adjust.method)
        result <- data.frame(gene = rownames(df),
            P.Value = P.Value, adj.P.Val = adj.P.Val)
    }
    return(result)
}


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
#' @param x matrix of Microarray data, each column is a sample, 
#' and each row is a gene. 
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
        x <- impute::impute.knn(x)$data
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
    x <- repAssign(x, string)
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
                rownames(clusterDf) <- 
                    paste(gene, seq(nrow(clusterDf)), sep = "_")
            } else {
                rownames(clusterDf) <- gene
            }
            
            colnames(clusterDf) <- colnames(mat)    
            matlist[[i]] <- clusterDf
        }
    }
    matlist <- do.call("rbind", matlist)
}