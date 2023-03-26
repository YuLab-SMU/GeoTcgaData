#' Do difference analysis of RNA-seq data
#'
#' @title differential_RNA
#' @rdname differential_RNA
#' @param counts a dataframe or numeric matrix of raw counts data, 
#' or SummarizedExperiment object
#' @param group sample groups
#' @param groupCol group column
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
#' @importFrom plyr rename
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment colData
#' @importFrom utils methods
#' @import methods
#' @import cqn
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
#' @importFrom plyr rename
#' @import cqn
differential_RNA <- function(counts, group, groupCol, method = "limma", 
                    geneLength = NULL,
                    gccontent = NULL, filter = TRUE, edgeRNorm = TRUE,
                    adjust.method = "BH", useTopconfects = TRUE, 
                    ucscData = FALSE) {
    method <- match.arg(method, c("DESeq2", "edgeR", "limma",
        "dearseq", "Wilcoxon", "NOISeq", "auto"))

    if (ucscData) {
        class(counts) <- "data.frame"
        counts[, 1] <- gsub("\\..*", "", counts[, 1])
        rownames(counts) <- counts[, 1]
        counts <- counts[, -1]
        counts <- round(2^counts) - 1
    }


    if (inherits(counts,  "SummarizedExperiment")) {
        se <- counts
        counts <- assays(se)$counts
        group <- colData(se)[, groupCol]
        names(group) <- rownames(colData(se))
    }


    cols <- !duplicated(colnames(counts))
    counts <- counts[, cols]
    group <- group[cols]
    if (min(table(group)) > 4 && method == "auto") {
        method <- "Wilcoxon"
    } else {
        method <- "limma"
    }

    ## use cqn to correct the bias
    correct<- TRUE
    uCovar <- NULL
    if (is.null(geneLength) || is.null(gccontent)) {
        correct<- FALSE
    } else {
        genes_gc <- intersect(names(geneLength), names(gccontent))
        uCovar <- data.frame(length = geneLength[genes_gc],
            gccontent = gccontent[genes_gc])
        rownames(uCovar) <- genes_gc
        counts <- counts[genes_gc, ]
    }
    d.mont <- edgeR::DGEList(counts = counts, group = group, genes = uCovar)
    if (filter) {
        keep <- edgeR::filterByExpr(d.mont)
        d.mont <- d.mont[keep, keep.lib.sizes = FALSE]
        counts <- counts[keep, ]
    }
    if (correct) {
        geneLength <- geneLength[rownames(counts)]
        gccontent <- gccontent[rownames(counts)]
        cqn.subset <- cqn::cqn(counts, lengths = geneLength, x = gccontent)
    }


    if (method == "DESeq2") {
        coldata <- data.frame(group)

        dds <- DESeq2::DESeqDataSetFromMatrix(
            countData = counts,
            colData = coldata, design = ~group
        )
        

        if (correct) {
            cqnOffset <- cqn.subset$glm.offset
            cqnNormFactors <- exp(cqnOffset)
            ## divide out the geometric mean
            ## https://support.bioconductor.org/p/89239/
            ## https://support.bioconductor.org/p/95683/
            normFactors <- cqnNormFactors / exp(rowMeans(log(cqnNormFactors)))
            DESeq2::normalizationFactors(dds) <- normFactors
        }
        DEGAll <- DESeq2::DESeq(dds)
        DEGAll_table <- NULL
        if (useTopconfects) {
            DEGAll_table <-
                topconfects::deseq2_confects(DEGAll, step = 0.05)$table
            rownames(DEGAll_table) <- DEGAll_table$name
        }
        DEGAll <- DEGAll |>
            DESeq2::results(pAdjustMethod = adjust.method) |>
            as.data.frame() |>
            rename(c("log2FoldChange" = "logFC")) |>
            rename(c("pvalue" = "P.Value")) |>
            rename(c("padj" = "adj.P.Val"))
        DEGAll$length <- geneLength[rownames(DEGAll)]
        DEGAll$gccontent <- gccontent[rownames(DEGAll)]
        # DEGAll <- DEGAll[, c(ncol(DEGAll) - 1,
        #     ncol(DEGAll), 1:(ncol(DEGAll) - 2))]
        DEGAll <- DEGAll[, c(ncol(DEGAll) - 1, ncol(DEGAll),
            seq_len(ncol(DEGAll) - 2))]
        if (!is.null(DEGAll_table)) {
            genes <- intersect(rownames(DEGAll), rownames(DEGAll_table))
            DEGAll <- cbind(DEGAll[genes, ], DEGAll_table[genes, ])
            DEGAll <- DEGAll[order(DEGAll$P.Value), ]
        }
    } else {
        if (correct) {
            ## with cqn, there is no need to normalize
            ## using the normalization tools
            ## from edgeR, such as calcNormFactors.
            d.mont$offset <- cqn.subset$glm.offset
        } else {
            ## TMM Normalization
            d.mont <- edgeR::calcNormFactors(d.mont, method = "TMM")
        }
        if (method == "edgeR") {
            # design <- stats::model.matrix(~ group)
            design <- stats::model.matrix(~ d.mont$sample$group)
            if (min(table(d.mont$sample$group)) > 1) {
                d.mont <- edgeR::estimateDisp(d.mont, design) |>
                    edgeR::estimateGLMCommonDisp(design = design)
                DEGAll <- edgeR::glmQLFit(d.mont, design = design)
                DEGAll_table <- NULL
                if (useTopconfects) {
                    DEGAll_table <- topconfects::edger_confects(DEGAll,
                        fdr = 0.05,
                        coef = ncol(DEGAll$design),
                        step = 0.05
                    )$table
                    rownames(DEGAll_table) <- DEGAll_table$name
                }
                # edgeR::topTags(n = nrow(d.mont$counts)) |>
                DEGAll <- DEGAll |>
                    edgeR::glmQLFTest(coef = ncol(DEGAll$design)) |>
                    edgeR::topTags(n = Inf, adjust.method = adjust.method) |>
                    as.data.frame() |>
                    rename(c("FDR" = "adj.P.Val")) |>
                    rename(c("PValue" = "P.Value"))
                if (!is.null(DEGAll_table)) {
                    genes <- intersect(rownames(DEGAll), rownames(DEGAll_table))
                    DEGAll <- cbind(DEGAll[genes, ], DEGAll_table[genes, ])
                    DEGAll <- DEGAll[order(DEGAll$P.Value), ]
                }
            } else {
                DEGAll <- edgeR::glmFit(d.mont, dispersion = 0)
                DEGAll <- DEGAll |>
                    edgeR::glmLRT(coef = ncol(DEGAll$design)) |>
                    edgeR::topTags(n = Inf, adjust.method = adjust.method) |>
                    as.data.frame() |>
                    rename(c("FDR" = "adj.P.Val")) |>
                    rename(c("PValue" = "P.Value"))
            }
        }

        if (method == "limma") {
            comparison <- paste(unique(group), collapse = "-")
            group <- factor(group)
            design <- stats::model.matrix(~ 0 + group)
            colnames(design) <- levels(group)
            contrast.matrix <- limma::makeContrasts(
                contrasts = comparison,
                levels = design
            )
            DEGAll <- limma::voom(d.mont, design = design, plot = FALSE) |>
                limma::lmFit(design)
            DEGAll_table <- NULL
            if (useTopconfects) {
                DEGAll_table <- topconfects::limma_confects(DEGAll,
                    coef = 1,
                    fdr = 0.05
                )$table
                rownames(DEGAll_table) <- DEGAll_table$name
            }
            DEGAll <- DEGAll |>
                limma::contrasts.fit(contrast.matrix) |>
                limma::eBayes() |>
                limma::topTable(number = Inf, adjust.method = adjust.method)
            if (!is.null(DEGAll_table)) {
                genes <- intersect(rownames(DEGAll), rownames(DEGAll_table))
                DEGAll <- cbind(DEGAll[genes, ], DEGAll_table[genes, ])
                DEGAll <- DEGAll[order(DEGAll$P.Value), ]
            }
        }

        if (method == "dearseq") {
            group[group == unique(group)[1]] <- 1
            group[group == unique(group)[2]] <- 2
            conditions <- matrix(as.numeric(group), ncol = 1)
            dearseqTest <- "asymptotic"
            if (edgeRNorm) {
                count_norm <- edgeR::cpm(d.mont, log = TRUE)
                DEGAll <- dearseq::dear_seq(
                    exprmat = count_norm, variables2test = conditions,
                    which_test = dearseqTest, parallel_comp = FALSE,
                    preprocessed = TRUE
                )
            } else {
                DEGAll <- dearseq::dear_seq(
                    exprmat = as.matrix(counts), variables2test = conditions,
                    which_test = dearseqTest, parallel_comp = FALSE,
                    preprocessed = FALSE,
                    padjust_methods = adjust.method
                )
            }
            DEGAll <- DEGAll$pvals |>
                rename(c("adjPval" = "adj.P.Val")) |>
                rename(c("rawPval" = "P.Value"))
        }

        if (method == "Wilcoxon") {
            count_norm <- edgeR::cpm(d.mont, log = TRUE) |> as.data.frame()
            pvalues <- rep(0, nrow(count_norm))

            count_disease <- as.matrix(count_norm[, group == unique(group)[1]])
            count_normal <- as.matrix(count_norm[, group == unique(group)[2]])
            for (i in seq_len(length(pvalues))) {
                pvalues[i] <- stats::wilcox.test(count_disease[i, ],
                    count_normal[i, ])$p.value
            }
            fdr <- stats::p.adjust(pvalues, method = adjust.method)
            DEGAll <- data.frame(P.Value = pvalues, adj.P.Val = fdr)
            rownames(DEGAll) <- rownames(count_norm)
        }

        if (method == "NOISeq") {
            conditions <- factor(group)
            data <- NOISeq::readData(data = counts,
                factors = as.data.frame(conditions))
            res <- NOISeq::noiseqbio(data,
                k = 0.5, norm = "tmm", factor = "conditions",
                random.seed = 12345, filter = 1, cv.cutoff = 100, cpm = 1
            )
            DEGAll <- NOISeq::degenes(res, q = 0, M = NULL) |>
                rename(c("prob" = "P.Value"))
            DEGAll$adj.P.Val <- DEGAll$P.Value
        }
    }
    if ("P.Value" %in% colnames(DEGAll)) {
        DEGAll <- DEGAll[!is.na(DEGAll[, "P.Value"]), ]
    }

    return(DEGAll)
}
