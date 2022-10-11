
#' a data.frame of gene expression data
#'
#' the rowname is gene symbols
#'
#' the columns are gene expression values
#'
#' @format A data.frame with 10779 rows and 2 column
#'
"geneExpress"

# ' a matrix for Converting gene symbol to entrez_id or ensembl_gene_id
# '
# ' the columns represent "symbol", "locus_group", "locus_type", 
# ' "entrez_id" and "ensembl_gene_id"
# '
# ' @format A matrix with 37647 rows and 5 column
# '
"hgnc"

#' a matrix of gene expression data in TCGA
#'
#' the first column represents the gene symbol
#'
#' the other columns represent the expression(FPKM) of genes
#'
#' @format A matrix with 10 rows and 10 column
#'
"profile"

#' a matrix of gene expression data in GEO
#'
#' the first column represents the gene symbol
#'
#' the other columns represent the expression of genes
#'
#' @format A matrix with 32 rows and 20 column
#'
"ventricle"

#' a matrix of gene expression data in TCGA
#'
#' the first column represents the gene symbol
#'
#' the other columns represent the expression(count) of genes
#'
#' @format A matrix with 100 rows and 150 column
#'
"kegg_liver"

#' a matrix of gene expression data in GEO
#'
#' the first column represents the gene symbol
#'
#' the other columns represent the expression of genes
#'
#' @format A matrix with 999 rows and 3 column
#'
"GSE66705_sample2"

#' a matrix of module name, gene symbols, and the number of gene symbols
#'
#' @format A matrix with 176 rows and 3 column
#'
"module"


# ' a matrix for Converting gene symbol.
# '
# '
# ' @format A matrix with 43547 rows and 52 column
# '
"hgnc_file"


#' a data.frame of gene length and GC content
#'
#'
#' @format A data.frame with 27341 rows and 2 column
#'
"gene_cov"
