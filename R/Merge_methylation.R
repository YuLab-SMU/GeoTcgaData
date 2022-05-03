#' Merge methylation data downloaded from TCGA
#'
#' @param dirr a string for the directory of methylation data download from tcga
#' useing the tools gdc
#' @return a matrix, a combined methylation expression spectrum matrix
#' @export
#'
#' @examples
#' merge_result <- Merge_methy_tcga(system.file(file.path("extdata","methy"),package="GeoTcgaData"))
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
    message("file",1," is over")
    for(i in 2:length(tcga_dir)) {
        message("file",i," is over")
		filePath <- file.path(dirr,tcga_dir[i])
        methyFile <- get_methy_df(filePath)
        methyResult[, i] <- methyFile[, 2]
        samples[i] <- colnames(methyFile)[2]
		gc()
    }
    colnames(methyResult) <- samples
    cpg_info <- methyFile[, -2]
    return(list(methyResult, cpg_info))
	# if (level == "cpg") return(jieguo)
    # file1 <- jieguo_old
    # genes <- as.matrix(methyFile[,6])
    # for(i in seq_len(dim(methyFile)[1])) {
    #     aa <- unlist(strsplit(genes[i],";"))
    #     bb <- unique(aa)
    #     genes[i] <- paste(bb, collapse=";")
    # }

    # file1[, 1] <- genes
	# file1 <- as.matrix(file1)
    # colnames(file1) <- file1[1, ]
	# file1 <- file1[-1, ]
	# rep1_result <- rep1(file1, ";")
    # ave_result <- gene_ave(rep1_result, k = 1)
}

#' Read methylated data file and turn it into data frame
#'
#' @param filePath Path of files
#' @noRd
get_methy_df <- function(filePath) {
    methyDir <- dir(filePath)
    for(j in 1:length(methyDir)) {
        if(length(grep("jhu-usc", methyDir[j]))>0) {
    		file_name <- file.path(filePath, dir(filePath)[j])
            sample <- unlist(strsplit(dir(filePath)[j],"\\."))[6]
        }
    }
    methyFile <- data.table::fread(file_name,header=T)
    class(methyFile) <- "data.frame"
    colnames(methyFile)[2] <- sample
    return(methyFile)
}

#' Get methylation difference gene
#'
#' @param cpgData data.frame of cpg beta value
#' @param sampleGroup vector of sample group
#' @param combineMethod method to combine the cpg pvalues
#' @export
methyDiff <- function(cpgData, sampleGroup, combineMethod = RobustRankAggreg::rhoScores) {
    if (class(cpgData) == "list") {
        cpgData <- cpgData[[1]]
    }
    cpgData <- as.matrix(cpgData)
    # Use KNN to fill in missing values
    data.m <- quiet(impute::impute.knn(cpgData)$data)
    # normalize data
    myNorm <- ChAMP::champ.norm(beta=data.m, rgSet = NULL, mset = NULL)
    # Identify Differential Methylation Positions (DMP)
    myDMP <- ChAMP::champ.DMP(beta = myNorm, pheno = sampleGroup, adjPVal = 1)
    myDMP <- as.data.frame(myDMP)
    # gene level difference analysis
    pvalues <- myDMP[, c(14, 4)]
    pvalues <- pvalues[pvalues[, 1] != "", ]
    gene_pvalue <- stats::aggregate(pvalues[, 2], by = list(pvalues[, 1]),
        FUN = combineMethod)
    # get logFC of genes
    myNorm2 <- myNorm[rownames(myDMP), ]
    myNorm2 <- stats::aggregate(myNorm2, by = list(myDMP[, 14]), FUN = mean)
    myNorm2 <- myNorm2[myNorm2[, 1] != "", ]
    rownames(myNorm2) <- myNorm2[, 1]
    myNorm2 <- myNorm2[, -1]
    groups <- sort(unique(sampleGroup))
    logFC <- rowMeans(myNorm2[, sampleGroup == groups[1]]) - rowMeans(myNorm2[, sampleGroup == groups[2]])
    gene_pvalue$logFC <- logFC[gene_pvalue[, 1]]
    colnames(gene_pvalue) <- c("gene", "pvalue", "logFC")
    return(gene_pvalue)
}

# from hadley wickham in "https://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html"
#' Suppressing output
#'
#' @param x some code
#' @noRd
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}


#' Title
#'
#' @param methy data.frame of the methylation data, which can be downloaded from UCSC Xena.
#' @param cpg_gene data.frame, the first coloumn is cpg id, the second coloumn is gene id.
#' @param group a vector of "0" and "1" for group of samples. If null, the samples were divided into two groups: disease and normal.
#' @param missing_value Method to  impute missing expression data, one of "zero" and "knn".
#' @importFrom dplyr `%>%`
#' @export
#'
#' @examples
#' \dontrun{
#' methy_file <- "TCGA.THCA.sampleMap_HumanMethylation450.gz"
#' methy <- data.table::fread(methy_file, sep = "\t", header = T)
#' library(ChAMP)
#' myImport <- champ.import(directory=system.file("extdata",package="ChAMPdata"))
#' myfilter <- champ.filter(beta=myImport$beta,pd=myImport$pd,
#'     detP=myImport$detP,beadcount=myImport$beadcount)
#' cpg_gene <- hm450.manifest.hg19[, c("probeID", "gene_HGNC")]
#' methyDiff_ucsc <- methyDiff_df(methy, cpg_gene)
#' }
methyDiff_ucsc <- function(methy, cpg_gene, group = NULL, missing_value = "knn") {
    class(methy) <- "data.frame"
    rownames(methy) <- methy[, 1]
    cpgs <- rownames(methy)
    methy <- methy[, -1]
    # dada <- substr(colnames(methy), 1, 12)
    # gaga <- table(dada)
    # co_group <- names(gaga[gaga > 1])
    # methy <- methy[, which(dada %in% co_group)]
 
    if (is.null(group)) {
        group <- lapply(colnames(methy), function(x) {
            strsplit(x, "-")[[1]][4]}) %>% unlist()
    
        group <- substring(group, 1,1)
    } 
    if (missing_value == "zero") {
        methy[is.na(methy)] <- 0
    } else {
        methy <- quiet(impute::impute.knn(as.matrix(methy))$data)
    }
    # rownames(methy) <- cpgs
    # methy_normal <- methy[, group == "1"]
    # methy_disease <- methy[, group == "0"]
    # # colnames(methy_disease) <- sapply(colnames(methy_disease), function(x) {
    # #     y <- strsplit(x, "-")[[1]][1:3]
    # #     paste(y, collapse = "-")
    # # })
    # # colnames(methy_disease) <- gsub("TCGA-", "", colnames(methy_disease))

    
    #  # To avoid reporting errors
    # merge_result <- cbind(methy_disease, methy_normal)
    merge_result <- as.data.frame(methy)
    merge_result$gene <- cpg_gene[rownames(merge_result), 2]
    merge_result <- merge_result[, c(ncol(merge_result), 1:(ncol(merge_result)-1))]

    merge_result <- merge_result[!is.na(merge_result$gene), ]


    merge_result$gene <- as.character(merge_result$gene)
    merge_result2 <- rep1(merge_result, ";")

    merge_result3 <- gene_ave(merge_result2)

    ## id conversion
    # haha <- id_conversion_vector("symbol", "entrez_id", rownames(merge_result3))
    # haha2 <- as.matrix(haha)
    # rownames(haha2) <- haha2[, 1]
    # merge_result3 <- merge_result3[haha2[, 1], ]
    # merge_result3 <- as.matrix(merge_result3)
    # rownames(merge_result3) <- gsub(" ", "", haha2[, 2])
    # merge_result3 <- merge_result3[!duplicated(rownames(merge_result3)), ]
    # merge_result3 <- merge_result3[!is.na(rownames(merge_result3)), ]
    ## use limma to do differential expression analysis
    methyDiff_limma(merge_result3, group = group)
    # group <- factor(group)
    # design <- model.matrix(~0 + group)
    # colnames(design) <- levels(group)
    # contrast.matrix <- makeContrasts(normal - cancer,
    #                                  levels=design)

    # fit <- lmFit(merge_result3, design)
    # fit2 <- contrasts.fit(fit, contrast.matrix)
    # fit2 <- eBayes(fit2)
    # topTable(fit2,adjust='fdr',coef=1,number=Inf)
}
#' methyDiff_limma
#'
#' @param df data.frame of the methylation data, which can be downloaded from UCSC Xena.
#' @param group a vector, group of samples.
#' @export
methyDiff_limma <- function(df, group) {
    groups <- unique(group)
    # if group is a numberic vector(even for c("0", "1")), will get errors.
    group <- gsub(groups[1], "nromal", group)
    group <- gsub(groups[2], "disease", group)
    design <- stats::model.matrix(~0 + factor(group))
    colnames(design) <- levels(factor(group))
    contrast.matrix <- limma::makeContrasts(contrasts = paste(colnames(design)[2:1], 
            collapse = "-"), levels = colnames(design))

    fit <- limma::lmFit(df, design)
    fit <- limma::contrasts.fit(fit, contrast.matrix)
    fit <- limma::eBayes(fit)
    limma::topTable(fit, adjust='BH', number=Inf)  

    ## contrasts.fit is not necessory
    # groups <- unique(group)
    # group <- gsub(groups[1], "nromal", group)
    # group <- gsub(groups[2], "disease", group)
    # design <- stats::model.matrix(~factor(group))
 
    # fit2 <- lmFit(df, design)
    # fit2 <- eBayes(fit2)
    # topTable(fit2,coef=2, adjust='BH') 
}
    
