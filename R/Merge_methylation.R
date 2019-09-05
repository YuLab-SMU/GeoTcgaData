#' Title
#'
#' @param dirr a string for the directory of methylation data download from tcga
#' useing the tools gdc
#' @return a matrix
#' @export
#'
#' @examples
#' mearge_result <- Merge_methy_tcga("direc")
Merge_methy_tcga <- function(dirr) {
    if(dirr != "direc") {
        tcga_dir <- dir(dirr)
        #dirr_l <- paste(dirr,tcga_dir[1],sep="\\")
		dirr_l <- file.path(dirr, tcga_dir[1])
        aa <- dir(dirr_l)
        for(j in 1:length(aa)) {
            if(length(grep("jhu-usc",aa[j]))>0) {
                #file_name <- paste(dirr_l,dir(dirr_l)[j],sep="\\")
				file_name <- file.path(dirr_l, dir(dirr_l)[j])
                sample_l <- unlist(strsplit(dir(dirr_l)[j],"\\."))[6]
            }
        }
		
        file_l <- data.table::fread(file_name,header=F)
        cpg <- file_l[,1]
        sample_l <- unlist(strsplit(dir(dirr_l)[1],"\\."))[6]
        beta <- file_l[,2]
        beta[1] <- sample_l
        jieguo <- cbind(cpg, beta)

        for(i in 2:length(tcga_dir)) {
            #dirr_l <- paste(dirr,tcga_dir[i],sep="\\")
			dirr_l <- file.path(dirr,tcga_dir[i])
            aa <- dir(dirr_l)
            for(j in 1:length(aa)) {
                if(length(grep("jhu-usc",aa[j]))>0) {
                    #file_name <- paste(dirr_l,dir(dirr_l)[j],sep="\\")
					file_name <- file.path(dirr_l,dir(dirr_l)[j])
                    sample_l <- unlist(strsplit(dir(dirr_l)[j],"\\."))[6]
                }
            }

            file_l <- data.table::fread(file_name,header=F)
            cpg <- file_l[,1]
            beta <- file_l[,2]
            beta[1] <- sample_l
            jieguo <- cbind(jieguo, beta)
        }
		
		
    file1 <- jieguo
    file2 <- file_l
    genes <- as.matrix(file2[,6])
    for(i in seq_len(dim(file2)[1])) {
        aa <- unlist(strsplit(genes[i],";"))
        bb <- unique(aa)
        genes[i] <- paste(bb,collapse=";")
    }

    file1[,1] <- genes
	file1 <- as.matrix(file1)
    colnames(file1) <- file1[1,]
	file1 <- file1[-1,]
	rep1_result <- rep1(file1,";")
    #rep1(file1,"methy_gene.txt",";")


    #ventricle <- data.table::fread("methy_gene.txt",sep="\t",header=F)
    #gene_ave(ventricle
    ave_result <- gene_ave(rep1_result)
    } else {message("please give your directory of methylation data!")}

}
