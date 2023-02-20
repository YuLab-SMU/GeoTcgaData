#' Get methylation difference gene
#'
#' @title differential_methy
#' @rdname differential_methy
#' @param cpgData data.frame of cpg beta value, , or SummarizedExperiment object
#' @param sampleGroup vector of sample group
#' @param combineMethod method to combine the cpg pvalues, 
#' a function or one of "stouffer", "fisher" and "rhoScores".
#' @param missing_value Method to impute missing expression data,
#' one of "zero" and "knn".
#' @param cpg2gene data.frame to annotate cpg locus to gene
#' @param normMethod Method to do normalization: "PBC" or "BMIQ".
#' @param region region of genes, one of "Body", "TSS1500", "TSS200",
#' "3'UTR", "1stExon", "5'UTR", and "IGR". Only used when cpg2gene is NULL.
#' @param model if "cpg", step1: calculate difference cpgs;
#' step2: calculate difference genes.
#' if "gene", step1: calculate the methylation level of genes;
#' step2: calculate difference genes.
#' @param adjust.method character string specifying the method
#' used to adjust p-values for multiple testing.
#' See \link{p.adjust} for possible values.
#' @param ucscData Logical, whether the data comes from UCSC Xena.
#' @param adjPvalCutoff adjusted pvalue cutoff
#' @param ... additional parameters
#' @return data.frame
#' @export
#' @examples
#' \donttest{
#' # use TCGAbiolinks data
#' library(TCGAbiolinks)
#' query <- GDCquery(project = "TCGA-ACC",
#'     data.category = "DNA Methylation",
#'     data.type = "Methylation Beta Value",
#'     platform = "Illumina Human Methylation 450")
#' GDCdownload(query, method = "api", files.per.chunk = 5,
#'     directory = Your_Path)
#' merge_result <- Merge_methy_tcga(Your_Path_to_DNA_Methylation_data)
#' library(ChAMP) # To avoid reporting errors
#' differential_gene <- differential_methy(cpgData = merge_result,
#'     sampleGroup = sample(c("C","T"),
#'     ncol(merge_result[[1]]), replace = TRUE))
#' }
#' # use user defined data
#' library(ChAMP)
#' cpgData <- matrix(runif(2000), nrow = 200, ncol = 10)
#' rownames(cpgData) <- paste0("cpg", seq_len(200))
#' colnames(cpgData) <- paste0("sample", seq_len(10))
#' sampleGroup <- c(rep("group1", 5), rep("group2", 5))
#' names(sampleGroup) <- colnames(cpgData)
#' cpg2gene <- data.frame(cpg = rownames(cpgData), 
#'     gene = rep(paste0("gene", seq_len(20)), 10))
#' result <- differential_methy(cpgData, sampleGroup, 
#'     cpg2gene = cpg2gene, normMethod = NULL)
#' # use SummarizedExperiment object input
#' library(ChAMP)
#' cpgData <- matrix(runif(2000), nrow = 200, ncol = 10)
#' rownames(cpgData) <- paste0("cpg", seq_len(200))
#' colnames(cpgData) <- paste0("sample", seq_len(10))
#' sampleGroup <- c(rep("group1", 5), rep("group2", 5))
#' names(sampleGroup) <- colnames(cpgData)
#' cpg2gene <- data.frame(cpg = rownames(cpgData), 
#'     gene = rep(paste0("gene", seq_len(20)), 10))
#' colData <- S4Vectors::DataFrame(
#'     row.names = colnames(cpgData),
#'     group = sampleGroup
#' )
#' data <- SummarizedExperiment::SummarizedExperiment(
#'          assays=S4Vectors::SimpleList(counts=cpgData),
#'          colData = colData)
#' result <- differential_methy(cpgData = data, 
#'     groupCol = "group", normMethod = NULL, 
#'     cpg2gene = cpg2gene)  
setGeneric("differential_methy",
            function(cpgData,  ...)
                standardGeneric("differential_methy")
            )


#' Do difference analysis of RNA-seq data
#'
#' @title differential_RNA
#' @rdname differential_RNA
#' @param counts a dataframe or numeric matrix of raw counts data, 
#' or SummarizedExperiment object
#' @param group sample groups
#' @param method one of "DESeq2", "edgeR" , "limma", "dearseq",
#' "NOISeq", "Wilcoxon", and "auto".
#' @param geneLength a vector of gene length.
#' @param gccontent a vector of gene GC content.
#' @param filter if TRUE, use filterByExpr to filter genes.
#' @param edgeRNorm if TRUE, use edgeR to do normalization for dearseq method.
#' @param adjust.method character string specifying the method used to
#' adjust p-values for multiple testing.
#' See \link{p.adjust} for possible values.
#' @param useTopconfects if TRUE, use topconfects to provide a
#'    more biologically useful ranked gene list.
#' @param ucscData Logical, whether the data comes from UCSC Xena.
#' @param ... additional parameters
#' @return data.frame
#' @export
#'
#' @examples
#' \donttest{
# use `TCGAbiolinks` to download TCGA data
#' library(TCGAbiolinks)
#'
#' query <- GDCquery(
#'     project = "TCGA-ACC",
#'     data.category = "Transcriptome Profiling",
#'     data.type = "Gene Expression Quantification",
#'     workflow.type = "STAR - Counts"
#' )
#'
#' GDCdownload(query,
#'     method = "api", files.per.chunk = 3,
#'     directory = Your_Path
#' )
#'
#' dataRNA <- GDCprepare(
#'     query = query, directory = Your_Path,
#'     save = TRUE, save.filename = "dataRNA.RData"
#' )
#' ## get raw count matrix
#' dataPrep <- TCGAanalyze_Preprocessing(
#'     object = dataRNA,
#'     cor.cut = 0.6,
#'     datatype = "STAR - Counts"
#' )
#'
#' # Use `differential_RNA` to do difference analysis.
#' # We provide the data of human gene length and GC content in `gene_cov`.
#' group <- sample(c("grp1", "grp2"), ncol(dataPrep), replace = TRUE)
#' library(cqn) # To avoid reporting errors: there is no function "rq"
#' ## get gene length and GC content
#' library(org.Hs.eg.db)
#' genes_bitr <- bitr(rownames(gene_cov),
#'     fromType = "ENTREZID", toType = "ENSEMBL",
#'     OrgDb = org.Hs.eg.db, drop = TRUE
#' )
#' genes_bitr <- genes_bitr[!duplicated(genes_bitr[, 2]), ]
#' gene_cov2 <- gene_cov[genes_bitr$ENTREZID, ]
#' rownames(gene_cov2) <- genes_bitr$ENSEMBL
#' genes <- intersect(rownames(dataPrep), rownames(gene_cov2))
#' dataPrep <- dataPrep[genes, ]
#' geneLength <- gene_cov2(genes, "length")
#' gccontent <- gene_cov2(genes, "GC")
#' names(geneLength) <- names(gccontent) <- genes
#' ##    Difference analysis
#' DEGAll <- differential_RNA(
#'     counts = dataPrep, group = group,
#'     geneLength = geneLength, gccontent = gccontent
#' )
#' # Use `clusterProfiler` to do enrichment analytics:
#' diffGenes <- DEGAll$logFC
#' names(diffGenes) <- rownames(DEGAll)
#' diffGenes <- sort(diffGenes, decreasing = TRUE)
#' library(clusterProfiler)
#' library(enrichplot)
#' library(org.Hs.eg.db)
#' gsego <- gseGO(gene = diffGenes, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
#' dotplot(gsego)
#' }
#' # use user-defined data
#' df <- matrix(rnbinom(400, mu = 4, size = 10), 25, 16)
#' df <- as.data.frame(df)
#' rownames(df) <- paste0("gene", 1:25)
#' colnames(df) <- paste0("sample", 1:16)
#' group <- sample(c("group1", "group2"), 16, replace = TRUE)
#' result <- differential_RNA(counts = df, group = group,
#'     filte = FALSE, method = "Wilcoxon")
#' # use SummarizedExperiment object input
#' df <- matrix(rnbinom(400, mu = 4, size = 10), 25, 16)
#' rownames(df) <- paste0("gene", 1:25)
#' colnames(df) <- paste0("sample", 1:16)
#' group <- sample(c("group1", "group2"), 16, replace = TRUE)
#' 
#' nrows <- 200; ncols <- 20
#'  counts <- matrix(
#'    runif(nrows * ncols, 1, 1e4), nrows,
#'    dimnames = list(paste0("cg",1:200),paste0("S",1:20))
#' )
#' 
#' colData <- S4Vectors::DataFrame(
#'   row.names = paste0("sample", 1:16),
#'   group = group
#' )
#' data <- SummarizedExperiment::SummarizedExperiment(
#'          assays=S4Vectors::SimpleList(counts=df),
#'          colData = colData)
#' 
#' result <- differential_RNA(counts = data, groupCol = "group",
#'     filte = FALSE, method = "Wilcoxon") 
setGeneric("differential_RNA",
            function(counts,  ...)
                standardGeneric("differential_RNA")
            )



