#' Do difference analysis of SNP data
#'
#' @param snpDf data.frame of SNP data.
#' @param sampleGroup vector of sample group.
#' @param method Method of combining the pvalue of multiple snp in a gene.
#' @export
diff_SNP <- function(snpDf, sampleGroup, method = min) {
    snpDf[!is.na(snpDf)] <- "mutation"
    snpDf[is.na(snpDf)] <- "wild"
    sampleGroup <- sampleGroup[!is.na(sampleGroup)]
    type1 <- which(sampleGroup == names(table(sampleGroup))[1])
    type2 <- which(sampleGroup == names(table(sampleGroup))[2])
    pvalue <- rep(0, nrow(snpDf))
    estimate <- rep(0, nrow(snpDf))
    for(i in seq_len(nrow(snpDf))) {
        type1_freq <- table(as.character(snpDf[i, type1]))
        type2_freq <- table(as.character(snpDf[i, type2]))
        df <- data.frame(type1 = as.numeric(type1_freq[c("wild", "mutation")]),
                         type2 = as.numeric(type2_freq[c("wild", "mutation")]))
        df[is.na(df)] <- 0
        fish <- stats::fisher.test(df)
        pvalue[i] <- fish$p.value
        estimate[i] <- fish$estimate
    }
    names(pvalue) <- names(estimate) <- sub("_.*", "", rownames(snpDf))
    pvalue <- stats::aggregate(pvalue, by = list(names(pvalue)), FUN = method)
    estimate <- stats::aggregate(estimate, by = list(names(estimate)), FUN = mean)
    return(data.frame(gene = pvalue[,1],  pvalue = pvalue[, 2], estimate = estimate[, 2]))
}

#' Do difference analysis of SNP data downloaded from TCGAbiolinks
#'
#' @param snpData data.frame of SNP data downloaded from TCGAbiolinks
#' @param sampleType vector of sample group
#' @export
#'
#' @examples
#' \dontrun{
#' library(TCGAbiolinks)
#' query <- GDCquery(project = "TCGA-ACC",
#'                   data.category = "Simple Nucleotide Variation",
#'                   data.type = "Masked Somatic Mutation",
#'                   workflow.type = "MuSE Variant Aggregation and Masking")
#' 
#' GDCdownload(query, method = "api", files.per.chunk = 5, directory = Your_Path)
#' 
#' data_snp <- GDCprepare(query = query, 
#'                    save = TRUE, 
#'                    directory =  "Your_Path") 
#' samples <- unique(data_snp$Tumor_Sample_Barcode)
#' sampleType <- sample(c("A","B"), length(samples), replace = TRUE)
#' names(sampleType) <- samples
#' pvalue <- diff_SNP_tcga(snpData = data_snp, sampleType = sampleType)
#' }
diff_SNP_tcga <- function(snpData, sampleType) {
    Tumor_Sample_Barcode <- Variant_Classification <- NULL
    snpName <- paste(snpData$Hugo_Symbol, snpData$Start_Position, sep = "_")
    snpData <- snpData[, c("Variant_Classification", "Tumor_Sample_Barcode")]
    snpData$snp <- snpName
    snpData <- tidyr::spread(snpData, Tumor_Sample_Barcode, Variant_Classification)
    snpData <- as.data.frame(snpData)
    i <- match(colnames(snpData), names(sampleType))
    sampleType <- sampleType[i]
    rownames(snpData) <- snpData$snp
    snpData <- snpData[, -1]
    pvalue <- diff_SNP(snpDf = snpData, sampleGroup = sampleType)
    return(pvalue)
}
