#' Do difference analysis of RNA-seq data
#'
#' @param counts a dataframe or numeric matrix of raw counts data
#' @param group sample groups
#' @param method one of "DESeq2", "edgeR" , "limma", "dearseq" and "Wilcoxon".
#' @param geneLength a vector of gene length.
#' @param gccontent a vector of gene GC content.
#' @param filter if TRUE, use filterByExpr to filter genes.
#' @param edgeRNorm if TRUE, use edgeR to do normalization for dearseq method.
#' @param adjust.method character string specifying the method used to adjust p-values for multiple testing. 
#' See \link{p.adjust} for possible values.
#' @importFrom magrittr %>%
#' @importFrom plyr rename
#' @import cqn
#' @export
#'
#' @examples
#' \dontrun{
# use `TCGAbiolinks` to download TCGA data
#' library(TCGAbiolinks)
#'                 
#' query <- GDCquery(project = "TCGA-ACC",
#'                   data.category = "Transcriptome Profiling",
#'                   data.type = "Gene Expression Quantification", 
#'                   workflow.type = "STAR - Counts")
#'                   
#' GDCdownload(query, method = "api", files.per.chunk = 3, 
#'     directory = Your_Path)
#' 
#' dataRNA <- GDCprepare(query = query, directory = Your_Path,
#'                       save = TRUE, save.filename = "dataRNA.RData")
#' ## get raw count matrix                         
#' dataPrep <- TCGAanalyze_Preprocessing(object = dataRNA,
#'                                       cor.cut = 0.6,
#'                                       datatype = "STAR - Counts")
#' 
#' # Use `diff_RNA` to do difference analysis. 
#' # We provide the data of human gene length and GC content in `gene_cov`.
#' group <- sample(c("grp1", "grp2"), ncol(dataPrep), replace = TRUE)
#' library(cqn) # To avoid reporting errors: there is no function "rq"
#' ## get gene length and GC content
#' library(org.Hs.eg.db)
#' genes_bitr <- bitr(rownames(gene_cov), fromType = "ENTREZID", toType = "ENSEMBL", 
#'          OrgDb = org.Hs.eg.db, drop = TRUE)
#' genes_bitr <- genes_bitr[!duplicated(genes_bitr[,2]), ]
#' gene_cov2 <- gene_cov[genes_bitr$ENTREZID, ]
#' rownames(gene_cov2) <- genes_bitr$ENSEMBL
#' genes <- intersect(rownames(dataPrep), rownames(gene_cov2))
#' dataPrep <- dataPrep[genes, ]
#' geneLength <- gene_cov2(genes, "length")
#' gccontent <- gene_cov2(genes, "GC")
#' names(geneLength) <- names(gccontent) <- genes
#' ##  Difference analysis
#' DEGAll <- diff_RNA(counts = dataPrep, group = group, 
#'                    geneLength = geneLength, gccontent = gccontent)
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
diff_RNA <- function(counts, group, method='limma', geneLength = NULL, 
                     gccontent = NULL, filter = TRUE, edgeRNorm = TRUE,
                     adjust.method = "BH") {

    method <- match.arg(method, c("DESeq2", "edgeR", "limma", "dearseq", "Wilcoxon"))
    cols <- !duplicated(colnames(counts))
    counts <- counts[, cols]
    group <- group[cols]
    ## use cqn to correct the bias
    correst <- TRUE
    uCovar <- NULL
    if (is.null(geneLength) || is.null(gccontent)) {
        correst <- FALSE
    } else {     
        genes_gc <- intersect(names(geneLength), names(gccontent))
        uCovar <- data.frame(length = geneLength[genes_gc], gccontent = gccontent[genes_gc])      
        rownames(uCovar) <- genes_gc
        counts <- counts[genes_gc, ]
    }
    d.mont <- edgeR::DGEList(counts = counts, group = group, genes = uCovar)
    if (filter) {      
        keep <- edgeR::filterByExpr(d.mont)
        d.mont <- d.mont[keep, keep.lib.sizes=FALSE]
        counts <- counts[keep, ]   
    }
    geneLength <- geneLength[rownames(counts)]   
    gccontent <- gccontent[rownames(counts)]
    if(correst) {
        cqn.subset <- cqn::cqn(counts, lengths = geneLength, x = gccontent)
    }
    if (method == 'DESeq2') {
        coldata <- data.frame(group)
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
            colData = coldata, design = ~ group)
        if (correst) {
            cqnOffset <- cqn.subset$glm.offset
            cqnNormFactors <- exp(cqnOffset)
            ## divide out the geometric mean
            ## https://support.bioconductor.org/p/89239/
            ## https://support.bioconductor.org/p/95683/
            normFactors <- cqnNormFactors / exp(rowMeans(log(cqnNormFactors)))
            DESeq2::normalizationFactors(dds) <- normFactors         
        }
        DEGAll <- DESeq2::DESeq(dds) %>% DESeq2::results(pAdjustMethod = adjust.method) %>% 
            as.data.frame() %>% 
            rename(c("log2FoldChange" = "logFC")) %>% 
            rename(c("pvalue" = "P.Value")) %>% 
            rename(c("padj" = "adj.P.Val"))
        DEGAll$length <- geneLength[rownames(DEGAll)]
        DEGAll$gccontent <- gccontent[rownames(DEGAll)]
        DEGAll <- DEGAll[, c(ncol(DEGAll)-1, ncol(DEGAll), 1:(ncol(DEGAll)- 2))]
    } else {       

        if (correst) {
            ## with cqn, there is no need to normalize using the normalization tools
            ## from edgeR, such as calcNormFactors. 
            d.mont$offset <- cqn.subset$glm.offset
        } else {
            ## TMM Normalization
            d.mont <- edgeR::calcNormFactors(d.mont)
        }
        if (method == "edgeR") {
            # design <- stats::model.matrix(~ group) 
            design <- stats::model.matrix(~ d.mont$sample$group) 
            d.mont <- edgeR::estimateDisp(d.mont, design)         
            DEGAll <- edgeR::estimateGLMCommonDisp(d.mont, design = design) %>%
                edgeR::glmFit(design = design) %>%
                edgeR::glmLRT(coef = 2) %>%
                # edgeR::topTags(n = nrow(d.mont$counts)) %>%
                edgeR::topTags(n = Inf, adjust.method = adjust.method) %>%
                as.data.frame() %>% 
                rename(c("FDR" = "adj.P.Val")) %>% 
                rename(c("PValue" = "P.Value"))
        }

        if(method == "limma") {
            comparison <- paste(unique(group), collapse = "-")
            group <- factor(group)
            design <- stats::model.matrix(~0+group)
            colnames(design) <- levels(group)
            contrast.matrix <- limma::makeContrasts(contrasts=comparison,
                levels=design)
            DEGAll <- limma::voom(d.mont, design=design, plot = FALSE) %>%
                limma::lmFit(design) %>%
                limma::contrasts.fit(contrast.matrix) %>%
                limma::eBayes() %>%
                limma::topTable(n = Inf, adjust.method = adjust.method)
        }

        if (method == "dearseq") {
            conditions <- matrix(as.numeric(group), ncol=1)
            dearseqTest <- "asymptotic"
            if(edgeRNorm){
                count_norm <- edgeR::cpm(d.mont, log=TRUE)              
                DEGAll <- dearseq::dear_seq(exprmat=count_norm, variables2test=conditions,
                    which_test=dearseqTest, parallel_comp=F, preprocessed=TRUE)
            }else{
                DEGAll <- dearseq::dear_seq(exprmat=as.matrix(counts), variables2test=conditions, 
                    which_test=dearseqTest, parallel_comp=F, preprocessed=FALSE,
                    padjust_methods = adjust.method)
            }
            DEGAll <- DEGAll$pvals %>% 
                rename(c("adjPval" = "adj.P.Val")) %>% 
                rename(c("rawPval" = "P.Value"))
        }

        if (method == "Wilcoxon") {
            count_norm <- edgeR::cpm(d.mont, log=TRUE) %>% as.data.frame()
            pvalues <- rep(0, nrow(count_norm))
            
            count_disease <- as.matrix(count_norm[, group == unique(group)[1]])
            count_normal <- as.matrix(count_norm[, group == unique(group)[2]])
            for (i in 1:length(pvalues)) {
               pvalues[i] <- stats::wilcox.test(count_disease[i, ], count_normal[i, ])$p.value
            }
            fdr <- stats::p.adjust(pvalues, method = adjust.method)
            DEGAll <- data.frame(P.Value = pvalues, adj.P.Val = fdr)
        }

        if (method == "NOISeq") {
            conditions <- factor(group)
            data <- NOISeq::readData(data=counts, factors=as.data.frame(conditions))
            res <- NOISeq::noiseqbio(data, k=0.5, norm="tmm", factor="conditions",
                random.seed = 12345, filter = 1, cv.cutoff = 100, cpm = 1)
            DEGAll <- NOISeq::degenes(res, q=0, M=NULL) %>% 
                rename(c("prob" = "P.Value"))
            DEGAll$adj.P.Val <- DEGAll$P.Value
        }
    }
    DEGAll <- DEGAll[!is.na(DEGAll[, "P.Value"]), ]
    return(DEGAll)
}

#' Do difference analysis of RNA-seq data downloaded from ucsc
#'
#' @param ucscfile a dataframe or numeric matrix of ucsc RNA-seq data
#' @param ... additional parameters
#' @export
#'
#' @examples

#' \dontrun{
#' ucscfile <- data.table::fread("TCGA-BRCA.htseq_counts.tsv.gz")
#' group <- sample(c("grp1", "grp2"), ncol(ucscfile) - 1, replace = TRUE)
#' result <- diff_RNA_ucsc(ucscfile, group = group)
#' }
diff_RNA_ucsc <- function(ucscfile, ...) {
    class(ucsc) <- "data.frame"
    ucsc[, 1] <- gsub("\\..*", "", ucsc[, 1])
    rownames(ucsc) <- ucsc[, 1]
    ucsc <- ucsc[, -1]
    ucsc <- round(2^ucsc) - 1
    diff_RNA(ucsc, ...)
}
