# deal_meta <- function(metadata_file) {
#     file1<-as.matrix(utils::read.table(metadata_file,sep="\t",header=F))

#     file2<-file1[grep(pattern="entity_id",file1[,1]),]
#     file3<-file1[grep(pattern="entity_submitter_id",file1[,1]),]

#     file3<-file3[nchar(file3)==57]

#     file2<-gsub("      entity_id: ","",file2)
#     file2<-gsub(", ","",file2)
#     file3<-gsub("      entity_submitter_id: ","",file3)
#     file3<-gsub(", ","",file3)
# 	return(cbind(file2,file3))
# }

# ' Merge the copy number variation data downloaded from TCGA using gdc
# '
# ' @param dirr a string of direction, catalogue of copy number variation data
# ' @param metadatafile a metadata file download from TCGA
# '
# ' @return a matrix,each column is a sample, each row is a gene
# ' @export
# '
# ' @examples
# ' metadatafile_name <- "metadata.cart.2018-11-09.json"
# ' \dontrun{result2 <- ann_merge(dirr = system.file(file.path("extdata","cnv"),
# ' package="GeoTcgaData"),metadatafile=metadatafile_name)}
# ann_merge <- function(dirr, metadatafile) {

#     tcga_dir <- dir(dirr)
#     output <- vector("list", length = length(tcga_dir))
# 	for(filek in seq_len(length(tcga_dir))) {
#         dirr_l <- file.path(dirr, tcga_dir[filek])
#         files <- dir(dirr_l)


#         for(j in seq_len(length(files))) {
#             if(length(grep(".txt",files[j]))>0) {

#                 file_name <- file.path(dirr_l, dir(dirr_l)[j])
#             }
#         }

#         # Each chromosome is compared separately to speed up
#         aa <- data.table::fread(file_name,header=TRUE)
#         class(aa) <- "data.frame"
#         aalist <- split(aa, aa$Chromosome)
#         genePoslist <- split(genePos, genePos$chr)
#         chrs <- intersect(names(aalist), names(genePoslist))
#         aalist <- aalist[chrs]
#         genePoslist <- genePoslist[chrs]
#         nlength <- unlist(lapply(aalist, nrow))
#         geneslist <- vector("list", length = length(chrs))
#         aalist2 <- vector("list", length = length(chrs))
#         for(i in seq_len(length(chrs))) {
#             aalisti <- aalist[[i]]
#             genePoslisti <- genePoslist[[i]]
#             genes <- vector("list", length = nlength[i])
#             for(j in seq_len(nrow(aalisti))) {
#                 rm1 <- genePoslisti$end < aalisti[j, "Start"]
#                 rm2 <- genePoslisti$start > aalisti[j, "End"]
#                 genes[[j]] <- genePoslisti[!(rm1 | rm2), "gene"]
#             }
#             genes <- unlist(lapply(genes, paste, collapse = ","))
#             aalist2[[i]] <- cbind(aalisti, genes)
#         }

#         # bb<-cbind(aa,genes)
#         bb <- do.call(rbind, aalist2)
#         bb <- bb[which(bb[,"genes"]!=""), ]
#         bb <- bb[, c("genes", "Segment_Mean")]

#         dd <- strsplit(bb$genes, ",")
#         ddLength <- unlist(lapply(dd, length))
#         output[[filek]] <- data.frame(genes = unlist(dd), 
#           Segment_Mean = rep(bb[, "Segment_Mean"], times = ddLength))
#         names(output)[filek] <- aa$GDC_Aliquot[1]
#     }

#     # genes_union <- unique(unlist(lapply(output, `[[`, 1)))
#     genes_union <- unique(unlist(lapply(output, function(x) x$genes)))
#     result <- matrix(0,length(genes_union),length(output))
#     rownames(result) <- genes_union
#     colnames(result) <- names(output)

#     for(file in names(output))
#     {
#         result[output[[file]]$genes, file] <- output[[file]]$Segment_Mean
#     }

#     metadata<-deal_meta(metadatafile)
#     rownames(metadata)<-metadata[,1]
#     colnames(result)<-metadata[colnames(result), 2]
#     resultdf <- as.data.frame(matrix(as.numeric(result), nrow(result)))
#     colnames(resultdf) <- colnames(result)
#     rownames(resultdf) <- rownames(result)
#     resultdf <- 2^(1+resultdf)
#     return(resultdf)
# }
