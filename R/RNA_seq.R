#' @rdname differential_RNA
#' @exportMethod differential_RNA
setMethod("differential_RNA", signature(counts = "matrix"),
            function(counts,  ...) {
                differential_RNA.matrix(counts,  ...)
            })


#' @rdname differential_RNA
#' @exportMethod differential_RNA
setMethod("differential_RNA", signature(counts = "data.frame"),
            function(counts,  ...) {
                differential_RNA.matrix(counts,  ...)
            })


#' @rdname differential_RNA
#' @exportMethod differential_RNA
setMethod("differential_RNA", signature(counts = "SummarizedExperiment"),
            function(counts,  ...) {
                differential_RNA.SummarizedExperiment(counts,  ...)
            })





#' @rdname differential_RNA
#'
#' @importFrom plyr rename
#' @import cqn
differential_RNA.matrix <- function(counts, group, method = "limma", 
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

#' @rdname differential_RNA
#' @param se SummarizedExperiment object
#' @param groupCol group column
#' @importFrom plyr rename
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment colData
#' @importFrom metap sumlog
#' @importFrom utils methods
#' @import methods
#' @import cqn
differential_RNA.SummarizedExperiment <- function(se, groupCol, 
                    method = "limma", geneLength = NULL,
                    gccontent = NULL, filter = TRUE, edgeRNorm = TRUE,
                    adjust.method = "BH", useTopconfects = TRUE) {
    # counts, group geneLength = NULL,gccontent = NULL
    counts <- assays(se)$counts
    group <- colData(se)[, groupCol]
    names(group) <- rownames(colData(se))
    differential_RNA(counts = counts, group = group, method = method, 
        geneLength = geneLength, gccontent = gccontent, 
        filter = filter, edgeRNorm = edgeRNorm,
        adjust.method = adjust.method, 
        useTopconfects = useTopconfects)

}

