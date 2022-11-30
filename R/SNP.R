#' Do difference analysis of SNP data
#'
#' @param snpDf data.frame of SNP data.
#' @param sampleGroup vector of sample group.
#' @param combineMethod Method of combining the
#' pvalue of multiple snp in a gene.
#' @return data.frame
#' @export
#'
#' @examples
#' \donttest{
#' library(TCGAbiolinks)
#' query <- GDCquery(
#'     project = "TCGA-CHOL",
#'     data.category = "Simple Nucleotide Variation",
#'     access = "open",
#'     legacy = FALSE,
#'     data.type = "Masked Somatic Mutation",
#'     workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
#' )
#' GDCdownload(query)
#' data_snp <- GDCprepare(query)
#' samples <- unique(data_snp$Tumor_Sample_Barcode)
#' sampleGroup <- sample(c("A", "B"), length(samples), replace = TRUE)
#' names(sampleGroup) <- samples
#' pvalue <- diff_SNP_tcga(snpData = data_snp, sampleGroup = sampleGroup)
#' }
#' # use demo data
#' snpDf <- matrix(sample(c("mutation", NA), 100, replace = TRUE), 10, 10)
#' snpDf <- as.data.frame(snpDf)
#' sampleGroup <- sample(c("A", "B"), 10, replace = TRUE)
#' result <- diff_SNP(snpDf, sampleGroup)
diff_SNP <- function(snpDf, sampleGroup, combineMethod = min) {
    snpDf[!is.na(snpDf)] <- "mutation"
    snpDf[is.na(snpDf)] <- "wild"
    sampleGroup <- sampleGroup[!is.na(sampleGroup)]
    type1 <- which(sampleGroup == names(table(sampleGroup))[1])
    type2 <- which(sampleGroup == names(table(sampleGroup))[2])
    pvalue <- rep(0, nrow(snpDf))
    estimate <- rep(0, nrow(snpDf))
    for (i in seq_len(nrow(snpDf))) {
        type1_freq <- table(as.character(snpDf[i, type1]))
        type2_freq <- table(as.character(snpDf[i, type2]))
        df <- data.frame(
            type1 = as.numeric(type1_freq[c("wild", "mutation")]),
            type2 = as.numeric(type2_freq[c("wild", "mutation")])
        )
        df[is.na(df)] <- 0
        fish <- stats::fisher.test(df)
        pvalue[i] <- fish$p.value
        estimate[i] <- fish$estimate
    }
    names(pvalue) <- names(estimate) <- sub("_.*", "", rownames(snpDf))
    if (!is.null(combineMethod)) {
        pvalue <- stats::aggregate(pvalue, by = list(names(pvalue)),
            FUN = combineMethod)
        estimate <- stats::aggregate(estimate,
            by = list(names(estimate)), FUN = mean)
        return(data.frame(gene = pvalue[, 1], pvalue = pvalue[, 2],
            estimate = estimate[, 2]))
    } else {
        return(data.frame(pvalue = pvalue, estimate = estimate))
    }
}

#' combine pvalues of SNP difference analysis result
#'
#' @param snpResult data.frame of SNP difference analysis result.
#' @param snp2gene data frame of two column: snp and gene.
#' @param combineMethod Method of combining the
#' pvalue of multiple snp in a gene.
#' @return data.frame
#' @export
combine_pvalue <- function(snpResult, snp2gene, combineMethod = min) {
        pvalue <- snpResult$pvalue
        estimate <- snpResult$estimate
        genes <- snp2gene[, 2]
        names(genes) <- snp2gene[, 1]
        snps <- rownames(snpResult)
        names(pvalue) <- names(estimate) <- genes[snps]
        pvalue <- stats::aggregate(pvalue, by = list(names(pvalue)),
            FUN = combineMethod)
        estimate <- stats::aggregate(estimate, by = list(names(estimate)),
            FUN = mean)
        return(data.frame(gene = pvalue[, 1], pvalue = pvalue[, 2],
            estimate = estimate[, 2]))

}

#' Do difference analysis of SNP data downloaded from TCGAbiolinks
#'
#' @param snpData data.frame of SNP data downloaded from TCGAbiolinks
#' @param sampleGroup vector of sample group
#' @param combineMethod Method of combining the pvalue of
#' multiple snp in a gene.
#' @return data.frame
#' @export
#'
#' @examples
#' \donttest{
#' library(TCGAbiolinks)
#' query <- GDCquery(
#'     project = "TCGA-CHOL",
#'     data.category = "Simple Nucleotide Variation",
#'     access = "open",
#'     legacy = FALSE,
#'     data.type = "Masked Somatic Mutation",
#'     workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
#' )
#' GDCdownload(query)
#' data_snp <- GDCprepare(query)
#' samples <- unique(data_snp$Tumor_Sample_Barcode)
#' sampleGroup <- sample(c("A", "B"), length(samples), replace = TRUE)
#' names(sampleGroup) <- samples
#' pvalue <- diff_SNP_tcga(snpData = data_snp, sampleGroup = sampleGroup)
#' }
#' # use demo data
#' snpDf <- matrix(sample(c("mutation", NA), 100, replace = TRUE), 10, 10)
#' snpDf <- as.data.frame(snpDf)
#' sampleGroup <- sample(c("A", "B"), 10, replace = TRUE)
#' result <- diff_SNP(snpDf, sampleGroup)
diff_SNP_tcga <- function(snpData, sampleGroup, combineMethod = NULL) {
    Tumor_Sample_Barcode <- Variant_Classification <- NULL
    snpName <- paste(snpData$Hugo_Symbol, snpData$Start_Position, sep = "_")
    ## 提取有用的信息
    snpData <- snpData[, c("Hugo_Symbol", "Start_Position", "Chromosome",
            "Variant_Classification", "Tumor_Sample_Barcode",
            "Variant_Type", "dbSNP_RS", "Mutation_Status",
            # "MAX_AF",
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Match_Norm_Validation_Allele1", "Match_Norm_Validation_Allele2")]


    snpData <- snpData[, c("Variant_Classification", "Tumor_Sample_Barcode")]
    snpData$snp <- snpName
    snpData <- tidyr::spread(snpData, Tumor_Sample_Barcode,
        Variant_Classification)
    snpData <- as.data.frame(snpData)
    i <- match(colnames(snpData), names(sampleGroup))
    sampleGroup <- sampleGroup[i]
    rownames(snpData) <- snpData$snp
    snpData <- snpData[, -1]
    pvalue <- diff_SNP(snpDf = snpData, sampleGroup = sampleGroup,
        combineMethod = combineMethod)
    return(pvalue)
}

#' Do difference analysis of SNP data downloaded from GEO
#'
#' @param snpData data.frame of SNP data downloaded from GEO
#' @param sampleGroup vector of sample group
#' @param method one of "Chisquare", "fisher",
#' and "CATT"(Cochran-Armitage trend test)
#' @return data.frame
#' @export
#' @examples
#' \donttest{
#' file1 <- read.table("GSE66903_series_matrix.txt.gz",
#'     fill=TRUE, comment.char="!", header = TRUE)
#' rownames(file1) <- file1[, 1]
#' snpData <- file1[, -1]
#' sampleGroup <- sample(c("A", "B"), ncol(snpData ), replace = TRUE)
#' names(sampleGroup) <- colnames(snpData)
#' snpData <- SNP_QC(snpData)
#' sampleGroup <- sample(c("A", "B"), ncol(snpData ), replace = TRUE)
#' result1 <- diff_SNP_GEO(snpData = snpData,
#'     sampleGroup = sampleGroup, method = "Chisquare")
#' }
#' # use demo data
#' snpDf <- matrix(sample(c("AA", "Aa", "aa"), 100, replace = TRUE), 10, 10)
#' snpDf <- as.data.frame(snpDf)
#' sampleGroup <- sample(c("A", "B"), 10, replace = TRUE)
#' result <- diff_SNP_GEO(snpDf, sampleGroup, method = "fisher")
diff_SNP_GEO <- function(snpData, sampleGroup, method = "Chisquare") {
    snpDf <- as.matrix(snpData)
    sampleGroup <- sampleGroup[!is.na(sampleGroup)]
    type1 <- which(sampleGroup == names(table(sampleGroup))[1])
    type2 <- which(sampleGroup == names(table(sampleGroup))[2])
    pvalue <- rep(1, nrow(snpDf))
    estimate <- rep(0, nrow(snpDf))
    for (i in seq_len(nrow(snpDf))) {
        type1_freq <- table(snpDf[i, type1])
        type2_freq <- table(snpDf[i, type2])
        types <- unique(snpDf[i, ])
        df <- data.frame(
            type1_freq = as.numeric(type1_freq[types]),
            type2_freq = as.numeric(type2_freq[types])
        )
        df[is.na(df)] <- 0
        if (nrow(df) > 2) {
            if (method == "fisher") {
                fish <- stats::fisher.test(df)
                pvalue[i] <- fish$p.value
                if (nrow(df) == 2) {
                        estimate[i] <- fish$estimate
                }
            }

            if (method == "Chisquare") {
                    pvalue[i] <- stats::chisq.test(df)$p.value
            }

            if(method == "CATT") {
                    pvalue[i] <- CATT::CATT(table = t(df))$p.value
            }
        }

    }
    names(pvalue) <- names(estimate) <- rownames(snpDf)

    return(data.frame(pvalue = pvalue, estimate = estimate))
}


get_maf <- function(x) {
    x <- x[x != "NoCall"]
    freq <- strsplit(x, split = "") |> unlist() |> table()
    min(freq) / sum(freq)
}

get_hwe <- function(x) {
    x <- x[x != "NoCall"]
    aa <- table(x)
    table_x <- as.numeric(aa)
    names(table_x) <- names(aa)
    # table_x <- table_x[sort(names(table_x))]
    freq <- strsplit(x, split = "") |> unlist() |> table()
    # freq <- freq[sort(names(freq))]
    table_y <- rep(0, 3)


    names(table_y) <- names(table_x)
    freq1 <- freq[1]/ sum(freq)
    freq2 <- freq[2]/ sum(freq)
    sum_freq <- length(x)
    table_y[paste0(names(freq)[1], names(freq)[1])] <- freq1 * freq1 * sum_freq
    if (length(table_x) > 1) {
        table_y[paste0(names(freq)[1],
            names(freq)[2])] <- freq1 * freq2 * sum_freq * 2
        table_y[paste0(names(freq)[2],
            names(freq)[1])] <- freq1 * freq2 * sum_freq * 2
    }
    if (length(table_x) > 2) {
            table_y[paste0(names(freq)[2],
                names(freq)[2])] <- freq2 * freq2 * sum_freq
    }
    df <- data.frame(table_x, table_y[names(table(x))])
    stats::chisq.test(df)$p.value
}

#' Do difference analysis of SNP data downloaded from TCGAbiolinks
#'
#' @param snpData data.frame of SNP data downloaded from TCGAbiolinks
#' @param geon filters out all variants with missing call rates
#' exceeding the provided value (default 0.02) to be removed
#' @param mind filters out all samples with missing call rates exceeding
#' the provided value (default 0.02) to be removed
#' @param maf filters out all variants with minor allele frequency below
#' the provided threshold
#' @param hwe filters out all variants which have Hardy-Weinberg
#' equilibrium exact test p-value below the provided threshold
#' @param miss character of miss value
#' @return data.frame
#' @export
#' @examples
#' # use demo data
#' snpDf <- matrix(sample(c("AA", "Aa", "aa"), 100, replace = TRUE), 10, 10)
#' snpDf <- as.data.frame(snpDf)
#' sampleGroup <- sample(c("A", "B"), 10, replace = TRUE)
#' result <- SNP_QC(snpDf)
SNP_QC <- function(snpData, geon = 0.02, mind = 0.02, maf = 0.05,
                    hwe = 1e-6, miss = "NoCall") {
    snpData_mat <- as.matrix(snpData)
    ## filter by 0.2
    aa <- snpData_mat |> apply(MARGIN = 2, FUN = function(x) {
        table(x)[miss] / length(x)
        }
    )

    aa[is.na(aa)] <- 0
    snpData_mat <- snpData_mat[, aa < 0.2]

    bb <- snpData_mat |> apply(MARGIN = 1, FUN = function(x) {
        table(x)[miss] / length(x)
        }
    )

    bb[is.na(bb)] <- 0
    snpData_mat <- snpData_mat[bb < 0.2, ]

    ## filter by cutoff
    aa <- snpData_mat |> apply(MARGIN = 2, FUN = function(x)
        {table(x)[miss] / length(x)})

    aa[is.na(aa)] <- 0
    snpData_mat <- snpData_mat[, aa < mind]

    bb <- snpData_mat |> apply(MARGIN = 1, FUN = function(x)
        {table(x)[miss] / length(x)})

    bb[is.na(bb)] <- 0
    snpData_mat <- snpData_mat[bb < geon, ]

    ## maf
    MAF <- snpData_mat |> apply(MARGIN = 1, FUN = get_maf)
    snpData_mat <- snpData_mat[MAF > maf, ]
    HWE <- snpData_mat |> apply(MARGIN = 1, FUN = get_hwe)
    snpData_mat <- snpData_mat[HWE > hwe, ] |> as.data.frame()
    return(snpData_mat)
}


