#' @rdname differential_methy
#' @exportMethod differential_methy
setMethod("differential_methy", signature(cpgData = "matrix"),
            function(cpgData,  ...) {
                differential_methy.matrix(cpgData,  ...)
            })


#' @rdname differential_methy
#' @exportMethod differential_methy
setMethod("differential_methy", signature(cpgData = "data.frame"),
            function(cpgData,  ...) {
                differential_methy.matrix(cpgData,  ...)
            })


#' @rdname differential_methy
#' @exportMethod differential_methy
setMethod("differential_methy", signature(cpgData = "list"),
            function(cpgData,  ...) {
                differential_methy.matrix(cpgData[[1]],  ...)
            })


#' @rdname differential_methy
#' @exportMethod differential_methy
setMethod("differential_methy", signature(cpgData = "SummarizedExperiment"),
            function(cpgData,  ...) {
                differential_methy.SummarizedExperiment(cpgData,  ...)
            })



#' @rdname differential_methy
#' @importFrom stats p.adjust
#' @importFrom metap sumz
#' @importFrom metap sumlog
differential_methy.matrix  <- function(cpgData, sampleGroup,
                    # combineMethod = RobustRankAggreg::rhoScores,
                    combineMethod = "stouffer",
                    missing_value = "knn", 
                    cpg2gene = NULL,
                    normMethod = "PBC", 
                    region = "TSS1500",
                    model = "gene",
                    adjust.method = "BH",
                    adjPvalCutoff = 0.05,
                    ucscData = FALSE) {
    region <- match.arg(region, c("Body", "TSS1500", "TSS200",
        "3'UTR", "1stExon", "5'UTR", "IGR"))
    model <- match.arg(model, c("cpg", "gene"))

    if (ucscData) {
        class(methy) <- "data.frame"
        rownames(methy) <- methy[, 1]
        cpgs <- rownames(methy)
        methy <- methy[, -1]
        group <- sampleGroup
        if (is.null(group)) {
            group <- lapply(colnames(methy), function(x) {
                strsplit(x, "-")[[1]][4]
            }) |> unlist()
    
            group <- substring(group, 1, 1)
        }
    }



    cpgData <- as.matrix(cpgData)
    # Use KNN to fill in missing values
    if (missing_value == "zero") {
        cpgData[is.na(cpgData)] <- 0
        data.m <- cpgData
    } else {
        data.m <- impute::impute.knn(cpgData)$data
    }

    # normalize data
    myNorm <- data.m
    if (!is.null(normMethod)) {
        myNorm <- ChAMP::champ.norm(beta = data.m, rgSet = NULL, 
            mset = NULL, method = normMethod)
    }
    if (!is.null(cpg2gene)) {
        cpg_gene <- cpg2gene
    } else {
        cpg_gene <- get_cpg_annotation(region = region)
    }


    if (model == "gene") {
        cpg_gene <-  split(cpg_gene[, 2], cpg_gene[, 1])   
        genes <- unlist(lapply(cpg_gene, function(x) {paste(x,collapse = ";")}))
        cpg_gene <- data.frame(cpg = names(cpg_gene), gene = genes)
        rownames(cpg_gene) <- cpg_gene[, 1]
        myNorm <- as.data.frame(myNorm)
        myNorm$gene <- cpg_gene[rownames(myNorm), 2]
        # myNorm <- myNorm[, c(ncol(myNorm), 1:(ncol(myNorm) - 1))]
        myNorm <- myNorm[, c(ncol(myNorm), seq_len(ncol(myNorm) - 1))]
        myNorm <- myNorm[!is.na(myNorm$gene), ]


        myNorm$gene <- as.character(myNorm$gene)
        myNorm2 <- repAssign(myNorm, ";")
        myNorm3 <- gene_ave(myNorm2)

        ## use limma to do differential expression analysis
        gene_pvalue <- differential_limma(myNorm3, group = sampleGroup,
            adjust.method = adjust.method)
        gene_pvalue$gene <- rownames(gene_pvalue)
    } else {
        # Identify Differential Methylation Positions (DMP)
        myDMP <- ChAMP::champ.DMP(beta = myNorm,
            pheno = sampleGroup, adjPVal = 1)
        myDMP <- as.data.frame(myDMP)

        # use cpg_gene to annotate CpGs
        pvalues <- cpg_gene
        pvalues$pvalue <- myDMP[cpg_gene[, 1], 4]
        # rownames(pvalues) <- pvalues[, 1]
        pvalues <- pvalues[!is.na(pvalues$pvalue), ]
        
        if (is.function(combineMethod)) {
            gene_pvalue <- stats::aggregate(pvalues[, 4],
                by = list(pvalues[, 2]),
                # FUN = combine_pvalue, combineMethod = combineMethod
                FUN = combineMethod
            )
            colnames(gene_pvalue) <- c("gene", "pvalue")
        } else {
            aa <- pvalues$pvalue
            bb <- split(aa, pvalues$gene)
            gene_pvalue <- data.frame(gene = names(bb), 
                pvalue = unlist(lapply(bb, function(x) x[1])))
            if (combineMethod == "stouffer") {
                
                myBetas <- myNorm[pvalues$cpg, ]
                myBetas <- split(as.data.frame(myBetas), pvalues$gene)
                correl <- lapply(myBetas, function(x) cor(t(x)))
                weights <- lapply(correl, function(x) 1/apply(x^2,1,sum))
                
                for (i in seq_len(nrow(gene_pvalue))) {
                    if (length(bb[[i]]) > 1) {
                        gene_pvalue[i, 2] <- sumz(bb[[i]], weights[[i]])$p
                    }       
                }
            }

            if (combineMethod == "fisher") {
                for (i in seq_len(nrow(gene_pvalue))) {
                    if (length(bb[[i]]) > 1) {
                        gene_pvalue[i, 2] <- sumlog(bb[[i]])$p
                    }       
                }
            }
        }
        


        # get logFC of genes
        myNorm2 <- myNorm[pvalues[, 1], ]
        myNorm2 <- stats::aggregate(myNorm2,
            by = list(pvalues[, 2]), FUN = mean)

        myNorm2 <- myNorm2[myNorm2[, 1] != "", ]
        rownames(myNorm2) <- myNorm2[, 1]
        myNorm2 <- myNorm2[, -1]
        groups <- sort(unique(sampleGroup))
        mean1 <- rowMeans(myNorm2[, sampleGroup == groups[1]], na.rm = TRUE)
        mean2 <- rowMeans(myNorm2[, sampleGroup == groups[2]], na.rm = TRUE)
        logFC <- mean1 - mean2            

        gene_pvalue$logFC <- logFC[gene_pvalue[, 1]]
        colnames(gene_pvalue) <- c("gene", "P.Value", "logFC")
        gene_pvalue$gene <- as.character(gene_pvalue$gene)
        gene_pvalue$adj.P.Val <- p.adjust(gene_pvalue$P.Value,
            method = adjust.method)
        rownames(gene_pvalue) <- gene_pvalue$gene
    }
    gene_pvalue <- gene_pvalue[gene_pvalue$adj.P.Val < adjPvalCutoff, ]
    return(gene_pvalue)                        
}


#' @rdname differential_methy
#' @param groupCol group column
#' @importFrom stats p.adjust
#' @importFrom metap sumz
#' @importFrom metap sumlog
differential_methy.SummarizedExperiment <- function(cpgData, groupCol,
                    combineMethod = "stouffer",
                    missing_value = "knn", 
                    cpg2gene = NULL,
                    normMethod = "PBC",
                    region = "TSS1500",
                    model = "gene",
                    adjust.method = "BH") {
    cpgData2 <- assays(cpgData)$counts
    group <- colData(cpgData)[, groupCol]
    names(group) <- rownames(colData(cpgData))
    differential_methy.matrix(cpgData = cpgData2, sampleGroup = group,
        combineMethod = combineMethod,
        missing_value = missing_value, 
        cpg2gene = cpg2gene,
        normMethod = normMethod,
        region = region,
        model = model,
        adjust.method = "BH")
}

#' differential_limma
#'
#' @param df data.frame of the omic data
#' @param group a vector, group of samples.
#' @param adjust.method adjust.method.
#' @return data.frame
#' @export
#' @examples
#' df <- matrix(runif(200), 25, 8)
#' df <- as.data.frame(df)
#' rownames(df) <- paste0("gene", 1:25)
#' colnames(df) <- paste0("sample", 1:8)
#' group <- sample(c("group1", "group2"), 8, replace = TRUE)
#' result <- differential_limma(df = df, group = group)
differential_limma <- function(df, group, adjust.method = "BH") {
    groups <- unique(group)
    # if group is a numberic vector(even for c("0", "1")), will get errors.
    group <- gsub(groups[1], "nromal", group)
    group <- gsub(groups[2], "disease", group)
    design <- stats::model.matrix(~ 0 + factor(group))
    colnames(design) <- levels(factor(group))
    contrast.matrix <- limma::makeContrasts(
        contrasts = paste(colnames(design)[2:1],
        collapse = "-"
    ), levels = colnames(design))

    fit <- limma::lmFit(df, design)
    fit <- limma::contrasts.fit(fit, contrast.matrix)
    fit <- limma::eBayes(fit)
    limma::topTable(fit, adjust.method = adjust.method, number = Inf)
    ## or limma::topTable(fit, coef = 1, adjust='BH', number=Inf)
    ## contrasts.fit is not necessory
    # groups <- unique(group)
    # group <- gsub(groups[1], "nromal", group)
    # group <- gsub(groups[2], "disease", group)
    # design <- stats::model.matrix(~factor(group))

    # fit2 <- lmFit(df, design)
    # fit2 <- eBayes(fit2)
    # topTable(fit2,coef=2, adjust='BH', number=Inf)

    ## coef parameter is not necessory：
    # opTable(fit2, adjust='BH', number=Inf)
}

#' Merge methylation data downloaded from TCGA
#'
#' When the methylation data is downloaded from TCGA, 
#' each sample is saved in a folder, which contains the methylation value file 
#' and the descriptive file. This function can directly 
#' extract and consolidate all folders.
#' @param dirr a string for the directory of methylation data download from tcga
#' useing the tools gdc
#' @return a matrix, a combined methylation expression spectrum matrix
#' @export
#'
#' @examples
#' merge_result <- Merge_methy_tcga(system.file(file.path("extdata", "methy"),
#'     package = "GeoTcgaData"))
Merge_methy_tcga <- function(dirr = NULL) {
    options(warn = -1)
    # file_num=1
    if (is.null(dirr)) stop("please give your directory of methylation data!")
    tcga_dir <- dir(dirr)
    filePath <- file.path(dirr, tcga_dir[1])
    methyFile <- get_methy_df(filePath)
    methyResult <- matrix(0, nrow = nrow(methyFile), ncol = length(tcga_dir))
    rownames(methyResult) <- methyFile[, "Composite Element REF"]
    samples <- rep(0, length(tcga_dir))
    methyResult[, 1] <- methyFile[, 2]
    samples[1] <- colnames(methyFile)[2]
    message("file", 1, " is over")
    for (i in 2:length(tcga_dir)) {
        message("file", i, " is over")
        filePath <- file.path(dirr, tcga_dir[i])
        methyFile <- get_methy_df(filePath)
        methyResult[, i] <- methyFile[, 2]
        samples[i] <- colnames(methyFile)[2]
        gc()
    }
    colnames(methyResult) <- samples
    cpg_info <- methyFile[, -2]
    return(list(methyResult = methyResult, cpg_info = cpg_info))
}

#' Read methylated data file and turn it into data frame
#'
#' @param filePath Path of files
#' @return data.frame
#' @noRd
get_methy_df <- function(filePath) {
    methyDir <- dir(filePath)
    for (j in seq_len(length(methyDir))) {
        if (length(grep("jhu-usc", methyDir[j])) > 0) {
            file_name <- file.path(filePath, dir(filePath)[j])
            sample <- unlist(strsplit(dir(filePath)[j], "\\."))[6]
        }
    }
    methyFile <- data.table::fread(file_name, header = TRUE)
    class(methyFile) <- "data.frame"
    colnames(methyFile)[2] <- sample
    return(methyFile)
}


get_cpg_annotation <- function(region = "TSS1500") {
    ## library to avoid errors.
    # library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    ann <- minfi::getAnnotation(
                IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19)
    ann <- as.data.frame(ann)
    cpg_gene <- ann[, c("Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group")]
    # cpg_gene <- cpg_gene[grep(region, cpg_gene$UCSC_RefGene_Group), ]
    cpg_gene <- cpg_gene[cpg_gene[, 2] != "", ]
    genelist <- strsplit(cpg_gene[, 2], ";")
    regionlist <- strsplit(cpg_gene[, 3], ";")
    geneLength <- unlist(lapply(genelist, length))
    cpgs <- rep(cpg_gene[, 1], times = geneLength)
    cpg_gene2 <- data.frame(cpg = cpgs, gene = unlist(genelist), 
        region = unlist(regionlist))
    cpg_gene2 <- cpg_gene2[grep(region, cpg_gene2$region), ]
    return(unique(cpg_gene2))
}
