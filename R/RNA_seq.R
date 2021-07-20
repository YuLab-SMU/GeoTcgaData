#' Do difference analysis of RNA-seq data
#'
#' @param counts a dataframe or numeric matrix of raw counts data
#' @param group sample groups
#' @param method one of "DESeq2", "edgeR" and "limma".
#' @param geneLength a vector of gene length.
#' @param gccontent a vector of gene GC content.
#' @importFrom magrittr %>%
#' @importFrom plyr rename
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
#'                   workflow.type = "HTSeq - Counts")
#'                   
#' GDCdownload(query, method = "api", files.per.chunk = 3, 
#'     directory = Your_Path)
#' 
#' dataRNA <- GDCprepare(query = query, directory = Your_Path,
#'                       save = TRUE, save.filename = "dataRNA.RData")
#' ## get raw count matrix                         
#' dataPrep <- TCGAanalyze_Preprocessing(object = dataRNA,
#'                                       cor.cut = 0.6,
#'                                       datatype = "HTSeq - Counts")
#' 
#' # Use `diff_RNA` to do difference analysis
#' ## Random value is used as gene length and GC content.
#' geneLength <- sample(1000:2000, nrow(dataPrep), replace = TRUE)
#' names(geneLength) <- colnames(dataPrep)
#' gccontent <- runif(nrow(dataPrep))
#' names(gccontent) <- colnames(dataPrep)
#' ## Random value is used as sample group.
#' group <- sample(c("grp1", "grp2"), ncol(dataPrep), replace = TRUE)
#' library(cqn) # To avoid reporting errors: there is no function "rq"
#' DEGAll <- diff_RNA(counts = dataPrep, group = group, 
#'                    geneLength = geneLength, gccontent = gccontent)
#' }
diff_RNA <- function(counts, group, method='DESeq2', geneLength = NULL, gccontent = NULL) {
    method <- match.arg(method, c("DESeq2", "edgeR", "limma"))

    ## use cqn to correct the bias
    correst <- TRUE
    if (is.null(geneLength) || is.null(gccontent)) {
        correst <- FALSE
    } else {
        cqn.subset <- cqn::cqn(counts, lengths = geneLength, x = gccontent)
    }

    if (method == 'DESeq2') {
        coldata <- data.frame(group)
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
            colData = coldata, design = ~ group)
        if (correst) {
            cqnOffset <- cqn.subset$glm.offset
            cqnNormFactors <- exp(cqnOffset)
            DESeq2::normalizationFactors(dds) <- cqnNormFactors
        }
        DEGAll <- DESeq2::DESeq(dds) %>% DESeq2::results() %>% as.data.frame() %>% 
            rename(c("log2FoldChange" = "logFC")) %>% 
            rename(c("pvalue" = "P.Value")) %>% 
            rename(c("padj" = "adj.P.Val"))
    } else {
        uCovar <- data.frame(length = geneLength, gccontent = gccontent)
        rownames(uCovar) <- rownames(counts)
        d.mont <- edgeR::DGEList(counts = counts, group = group, genes = uCovar)
        if (method == "edgeR") {
            design <- stats::model.matrix(~ group)
            if (correst) {
                d.mont$offset <- cqn.subset$glm.offset
            }
            DEGAll <- edgeR::estimateGLMCommonDisp(d.mont, design = design) %>%
                edgeR::glmFit(design = design) %>%
                edgeR::glmLRT(coef = 2) %>%
                edgeR::topTags(n = nrow(d.mont$counts)) %>%
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
                limma::topTable(coef=1, n = Inf)
        }
    }
    return(DEGAll)
}
