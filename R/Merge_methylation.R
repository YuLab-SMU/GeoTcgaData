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
