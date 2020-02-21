

# GeoTcgaData

The goal of GeoTcgaData is to deal with RNA-seq, DNA Methylation, and Copy 
number variation data in GEO and TCGA.

## :writing_hand: Authors
Erqiang Hu

College of Bioinformatics Science and Technology, Harbin Medical University

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/GeoTcgaData?color=green)](https://cran.r-project.org/package=GeoTcgaData)
![](https://cranlogs.r-pkg.org/badges/grand-total/GeoTcgaData?color=green)
![](https://cranlogs.r-pkg.org/badges/GeoTcgaData?color=green)
![](https://cranlogs.r-pkg.org/badges/last-week/GeoTcgaData?color=green)

## :arrow\_double\_down: Installation

Get the development version from github:

```r
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("huerqiang/GeoTcgaData")
```
Or  the released version from CRAN:

``` r
install.packages("GeoTcgaData")
```
GEO and TCGA provide us with a wealth of data, such as RNA-seq, DNA Methylation,  and Copy number variation data. It's easy to download data from TCGA using the  gdc tool, but processing these data into a format suitable for bioinformatics  analysis requires more work. This R package was developed to handle these data.

## Example

This is a basic example which shows you how to solve a common problem:

## RNA-seq data integration and differential gene extraction
The function `classify_sample` and `diff_gene` could get the differentially expressioned genes using `DESeq2` package. For examples:
```r
library(DESeq2)
profile2 <- classify_sample(kegg_liver) 
jieguo <- diff_gene(profile2)
```

The parameter ` kegg_liver` is a matrix or data.frame of gene expression data(count) in TCGA.

## DNA Methylation data integration 
The function Merge_methy_tcga could Merge methylation data downloaded from TCGA. This makes it easier to extract differentially methylated genes in the downstream analysis. For example:

```r
dirr = system.file(file.path("extdata","methy"),package="GeoTcgaData")
merge_result <- Merge_methy_tcga(dirr)
```

## Copy number variation data integration and differential gene extraction
The function `ann_merge` could merge the copy number variation data downloaded from TCGA using gdc. For example:

```r
metadatafile_name <- "metadata.cart.2018-11-09.json"
jieguo2 <- ann_merge(dirr = system.file(file.path("extdata","cnv"),package="GeoTcgaData"),metadatafile=metadatafile_name)
```

The parameter `dirr` is a string for the direction of copy number variation data downloaded from TCGA. The parameter `metadatafile` is the metadata file download from TCGA.
The function `prepare_chi` and `differential_cnv` could do chi-square test to find copy number variation differential genes. For example:

```r
jieguo3 <- matrix(c(-1.09150,-1.47120,-0.87050,-0.50880,
-0.50880,2.0,2.0,2.0,2.0,2.0,2.601962,2.621332,2.621332,
                    2.621332,2.621332,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,
                    2.0,2.0,2.0,2.0,2.0,2.0,2.0),nrow=5)
rownames(jieguo3) <- c("AJAP1", "FHAD1", "CLCNKB", "CROCCP2", "AL137798.3")
colnames(jieguo3) <- c("TCGA-DD-A4NS-10A-01D-A30U-01", "TCGA-ED-A82E-01A-11D-A34Y-01", 
"TCGA-WQ-A9G7-01A-11D-A36W-01", "TCGA-DD-AADN-01A-11D-A40Q-01", 
"TCGA-ZS-A9CD-10A-01D-A36Z-01", "TCGA-DD-A1EB-11A-11D-A12Y-01")
 rt <- prepare_chi(jieguo3)
 chiResult <- differential_cnv(rt)
```

The parameter of `prepare_chi` is the result of function `ann_merge` and the parameter of function `differential_cnv` is the result of prepare_chi.

## GEO chip data processing
The function `gene_ave` could average the expression data of different ids for the same gene in the GEO chip data. For example:

```r
aa <- c("Gene Symbol", "MARCH1", "MARC1", "MARCH1", "MARCH1", "MARCH1")
bb <- c("GSM1629982", "2.969058399", "4.722410064", "8.165514853", "8.24243893", "8.60815086")
cc <- c("GSM1629982", "3.969058399", "5.722410064", "7.165514853", "6.24243893", "7.60815086")
file1 <- data.frame(aa=aa,bb=bb,cc=cc)
result <- gene_ave(file1)
```

Multiple genes symbols may correspond to a same chip id. The result of function `rep1` is to assign the expression of this id to each gene, and function `rep2` deletes the expression. For example:

```r
aa <- c("MARCH1 /// MMA","MARC1","MARCH2 /// MARCH3",
        "MARCH3 /// MARCH4","MARCH1")
bb <- c("2.969058399","4.722410064","8.165514853","8.24243893","8.60815086")
cc <- c("3.969058399","5.722410064","7.165514853","6.24243893","7.60815086")
input_fil <- data.frame(aa=aa,bb=bb,cc=cc)
rep1_result <- rep1(input_fil," /// ")
rep2_result <- rep2(input_fil," /// ")
```

## Other downstream analyses

1. The function `id_conversion_vector` could convert gene id from one of `symbol`, `RefSeq_ID`, `Ensembl_ID`, `NCBI_Gene_ID`, `UCSC_ID`, and `UniProt_ID` , etc. to another. Use `id_ava()` to get all the convertible ids. For example:

```r
id_conversion_vector("symbol", "ensembl_gene_id", c("A2ML1", "A2ML1-AS1", "A4GALT", "A12M1", "AAAS")) 

```


When the user converts the Ensembl ID to other ids, the version number needs to be removed. For example, "ENSG00000186092.4" doesn't work, you need to change it to "ENSG00000186092".

Especially, the function id_conversion could convert  ENSEMBL gene id to gene Symbol in TCGA. For example:

```r
result <- id_conversion(profile)
```

The parameter profile is a data.frame or matrix of gene expression data in TCGA.

2. The function `countToFpkm_matrix` and `countToTpm_matrix` could convert count data to FPKM or TPM data.

```r
lung_squ_count2 <- matrix(c(1,2,3,4,5,6,7,8,9),ncol=3)
rownames(lung_squ_count2) <- c("DISC1","TCOF1","SPPL3")
colnames(lung_squ_count2) <- c("sample1","sample2","sample3")
jieguo <- countToFpkm_matrix(lung_squ_count2)
```

```r
lung_squ_count2 <- matrix(c(11,22,23,14,15,6,17,18,29),ncol=3)
rownames(lung_squ_count2) <- c("DISC1","TCOF1","SPPL3")
colnames(lung_squ_count2) <- c("sample1","sample2","sample3")
jieguo <- countToTpm_matrix(lung_squ_count2)
```
3. The function `tcga_cli_deal` could combine clinical information obtained from TCGA and extract survival data. For example:

```r
tcga_cli <- tcga_cli_deal(system.file(file.path("extdata","tcga_cli"),package="GeoTcgaData"))
```
