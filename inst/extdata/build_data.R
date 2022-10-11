#setwd("E:\\GeoTcgaData_work")
hgnc_file <- data.table::fread("E:\\GeoTcgaData_work\\hgnc_complete_set.txt", 
  sep = "\t", header = TRUE)
hgnc_file <- dplyr::select(hgnc_file, -c("alias_symbol", "alias_name", 
  "prev_symbol", "lsdb", "agr"))
class(hgnc_file) <- "data.frame"
gene_loc_len <- GeoTcgaData:::gene_loc_len
hgnc <- GeoTcgaData:::hgnc
# genePos <- GeoTcgaData:::genePos
# genePos <- as.data.frame(genePos)
# genePos$start <- as.integer(genePos$start)
# genePos$end <- as.integer(genePos$end)
# genePos$gene_len <- as.integer(genePos$gene_len)
genePos <- GeoTcgaData:::genePos
hgnc_file <- GeoTcgaData:::hgnc_file
# usethis::use_data(hgnc_file, gene_loc_len, hgnc, genePos, internal = TRUE,
#   compress = "xz", overwrite = TRUE)
usethis::use_data(hgnc_file, hgnc, gene_loc_len, 
  internal = TRUE, compress = "xz", overwrite = TRUE)


## gene_cov
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
hg38_TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
hg38 <- BSgenome.Hsapiens.UCSC.hg38
calc_gene_cov <- function(TxDB, BSGENOME){
    Gene <- genes(TxDB, single.strand.genes.only = FALSE)
    Exon <- exons(x = TxDB)
    Overlap <- findOverlaps(Exon, Gene)
    Exon <- Exon[queryHits(Overlap)]
    mcols(Exon)$gene_id <- mcols(Gene[subjectHits(Overlap)])$gene_id
    Exon <- split(Exon, mcols(Exon)$gene_id)
    Exon <- reduce(Exon)
    calculate_cov <- function(x){
        xlen <- sum(width(x))
        xseq <- BSgenome::getSeq(BSGENOME, x)
        xGC <- sum(Biostrings::letterFrequency(xseq, 'GC'))/xlen
        c(xlen, xGC)
    }
    gene_cov <- lapply(Exon, calculate_cov)
    gene_cov <- gene_cov[names(Gene)]
    gene_cov <- t(as.data.frame(gene_cov))
    rownames(gene_cov) <- names(Gene)
    colnames(gene_cov) <- c('length', 'GC')
    as.data.frame(gene_cov)
}
gene_cov <- calc_gene_cov(TxDB = hg38_TxDb, BSGENOME = hg38)



