---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```
# GeoTcgaData

<!-- badges: start -->
<!-- badges: end -->

The goal of GeoTcgaData is to deal with RNA-seq, DNA Methylation, and Copy 
number variation data in GEO and TCGA.

## Installation

You can install the released version of GeoTcgaData from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("GeoTcgaData")
```
GEO and TCGA provide us with a wealth of data, such as RNA-seq, DNA Methylation, 
and Copy number variation data. It's easy to download data from TCGA using the 
gdc tool, but processing these data into a format suitable for bioinformatics 
analysis requires more work. This R package was developed to handle these data.
## Example

This is a basic example which shows you how to solve a common problem:

```{r example}
library(GeoTcgaData)

#' This function can find the mean value of the gene in each module.
#'
#' param geneExpress a data.frame
#' param module a data.frame
#' param result a string
#'
#' @return a matrix, means the mean of gene expression value in
#' the same module
result <- cal_mean_module(geneExpress,module)

#' Average the expression data of different ids for the same gene in the chip
#' expression profile of GEO or TCGA
#'
#' param file1 a data.frame
#' param k a number
#'
#' @return a data.frame, the values of same genes in gene expression profile

#'
#' examples
aa <- c("Gene Symbol","MARCH1","MARC1","MARCH1","MARCH1","MARCH1")
bb <- c("GSM1629982","2.969058399","4.722410064","8.165514853","8.24243893","8.60815086")
cc <- c("GSM1629982","3.969058399","5.722410064","7.165514853","6.24243893","7.60815086")
file1 <- data.frame(aa=aa,bb=bb,cc=cc)
result <- gene_ave(file1)

#' Get the differentially expressioned genes using DESeq2 package 
#'
#' param profile a data.frame
#'
#' @return a data.frame, a intermediate results of DESeq2
#'
#' examples
profile2 <- classify_sample(profile)


#' Get the differentially expressioned genes using DESeq2 package
#'
#' param profile2 a result of classify_sample
#'
#' @return a matrix, information of differential expression genes
#'
#' examples
profile2 <- classify_sample(profile)
jieguo <- diff_gene(profile2)

#' Merge methylation data downloaded from TCGA
#'
#' param dirr a string for the directory of methylation data download from tcga
#' useing the tools gdc
#' @return a matrix, a combined methylation expression spectrum matrix
#'
#' examples
mearge_result <- Merge_methy_tcga("direc")

#' Combine clinical information obtained from TCGA and extract survival data
#'
#' param Files_dir1 a dir data
#'
#' @return a matrix, survival time and survival state in TCGA
#'
#' examples
tcga_cli_deal("your_clinical_directory")

#' Multiple genes symbols may correspond to a same id. Some people think 
#' that the expression value of this id should be 
#' given to each gene, and some people think that the expression value of
#' this id should be deleted. The result of rep1 is to assign the expression
#' of this id to each gene, and rep2 deletes the expression.
#'
#' param file1 input file, a data.frame or a matrixg
#' param string a string,sep of the gene
#'
#' return a data.frame,rep1 is to assign the expression
#' of this id to each gene, and rep2 deletes the expression.
#'
#' examples
aa <- c("MARCH1 /// MMA","MARC1","MARCH2 /// MARCH3",
        "MARCH3 /// MARCH4","MARCH1")
bb <- c("2.969058399","4.722410064","8.165514853","8.24243893","8.60815086")
cc <- c("3.969058399","5.722410064","7.165514853","6.24243893","7.60815086")
input_fil <- data.frame(aa=aa,bb=bb,cc=cc)
rep1_result <- rep1(input_fil," /// ")
rep1_result <- rep2(input_fil," /// ")

#' Convert  ENSEMBL gene id to gene Symbol in TCGA
#'
#' param profile a data.frame
#'
#' @return a data.frame, gene symbols and their expression value
#'
#' examples
result <- id_conversion(profile)



```

```{r, eval=FALSE, message=FALSE, warning=FALSE}
#' Title gene id conversion
#'
#' param from one of "symbol","RefSeq_ID","Ensembl_ID","NCBI_Gene_ID","UCSC_ID","UniProt_ID"
#' param to one of "symbol","RefSeq_ID","Ensembl_ID","NCBI_Gene_ID","UCSC_ID","UniProt_ID"
#' param IDs the gene id which needed to convert
#'
#' return a vector of genes
#' export
#'
#' examples
id_conversion_vector("symbol","Ensembl_ID",c("A2ML1","A2ML1-AS1","A4GALT","A12M1","AAAS"))

```

```{r, eval=FALSE, message=FALSE, warning=FALSE}

#' Title Convert count to FPKM
#'
#' @param counts_matrix a matrix, colnames of counts_matrix are sample name,
#' rownames of counts_matrix are gene symbols
#'
#' @return  a matrix
#' @export
#'
#' @examples
#' lung_squ_count2 <- matrix(c(1,2,3,4,5,6,7,8,9),ncol=3)
#' rownames(lung_squ_count2) <- c("DISC1","TCOF1","SPPL3")
#' colnames(lung_squ_count2) <- c("sample1","sample2","sample3")
jieguo <- countToFpkm_matrix(lung_squ_count2)

```


```{r, eval=FALSE, message=FALSE, warning=FALSE}

#' Title Convert fpkm to Tpm
#'
#' @param fpkm_matrix a matrix, colnames of fpkm_matrix are sample name,
#' rownames of fpkm_matrix are gene symbols
#'
#' @return a matrix
#' @export
#'
#' @examples
#' lung_squ_count2 <- matrix(c(0.11,0.22,0.43,0.14,0.875,0.66,0.77,0.18,0.29),ncol=3)
#' rownames(lung_squ_count2) <- c("DISC1","TCOF1","SPPL3")
#' colnames(lung_squ_count2) <- c("sample1","sample2","sample3")
jieguo <- countToTpm_matrix(lung_squ_count2)

```



```{r, eval=FALSE, message=FALSE, warning=FALSE}

#' Title merge the copy number variation data downloaded from TCGA using gdc
#'
#' @param dirr a string of direction, catalogue of copy number variation data
#' @param metadatafile a metadata file download from TCGA
#'
#' @return a matrix,each column is a sample, each row is a gene
#' @export
#'
#' @examples
#' metadatafile_name <- "metadata.cart.2018-11-09.json"
jieguo2 <- ann_merge(dirr = system.file(file.path("extdata","cnv"),package="GeoTcgaData"),metadatafile=metadatafile_name)

```

```{r, eval=FALSE, message=FALSE, warning=FALSE}

#' Title preparer file for chi-square test
#'
#' @param jieguo2 result of ann_merge()
#'
#' @return a matrix
#' @export
#'
#' @examples
jieguo3 <- matrix(c(-1.09150,-1.47120,-0.87050,-0.50880,
                    -0.50880,2.0,2.0,2.0,2.0,2.0,2.601962,2.621332,2.621332,
                     2.621332,2.621332,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,
                     2.0,2.0,2.0,2.0,2.0,2.0,2.0),nrow=5)
rownames(jieguo3) <- c("AJAP1","FHAD1","CLCNKB","CROCCP2","AL137798.3")
colnames(jieguo3) <- c("TCGA-DD-A4NS-10A-01D-A30U-01","TCGA-ED-A82E-01A-11D-A34Y-01",
"TCGA-WQ-A9G7-01A-11D-A36W-01","TCGA-DD-AADN-01A-11D-A40Q-01",
"TCGA-ZS-A9CD-10A-01D-A36Z-01","TCGA-DD-A1EB-11A-11D-A12Y-01")
cnv_chi_file <- prepare_chi(jieguo3)

```


```{r, eval=FALSE, message=FALSE, warning=FALSE}

#' Title do chi-square test to find differential genes
#'
#' @param cnv_chi_file result of prepare_chi()
#'
#' @return a matrix
#' @export
#'
#' @examples
jieguo3 <- matrix(c(-1.09150,-1.47120,-0.87050,-0.50880,
                     -0.50880,2.0,2.0,2.0,2.0,2.0,2.601962,2.621332,2.621332,
                    2.621332,2.621332,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,
                    2.0,2.0,2.0,2.0,2.0,2.0,2.0),nrow=5)
 rownames(jieguo3) <- c("AJAP1","FHAD1","CLCNKB","CROCCP2","AL137798.3")
colnames(jieguo3) <- c("TCGA-DD-A4NS-10A-01D-A30U-01","TCGA-ED-A82E-01A-11D-A34Y-01",
"TCGA-WQ-A9G7-01A-11D-A36W-01","TCGA-DD-AADN-01A-11D-A40Q-01",
"TCGA-ZS-A9CD-10A-01D-A36Z-01","TCGA-DD-A1EB-11A-11D-A12Y-01")
 rt <- prepare_chi(jieguo3)
 chiResult <- differential_cnv(rt)

```
