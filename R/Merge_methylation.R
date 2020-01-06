#' Merge methylation data downloaded from TCGA
#'
#' @param dirr a string for the directory of methylation data download from tcga
#' useing the tools gdc
#' @return a matrix, a combined methylation expression spectrum matrix
#' @export
#'
#' @examples
#' merge_result <- Merge_methy_tcga(system.file(file.path("extdata","methy"),package="GeoTcgaData"))
Merge_methy_tcga <- function(dirr) {
    options(warn = -1)
    file_num=1
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
		
        file_l <- data.table::fread(file_name,header=FALSE)
		file_ll <- file_l
		
        cpg <- file_l[,1]
        sample_l <- unlist(strsplit(dir(dirr_l)[1],"\\."))[6]
        betaa <- file_l[,2]
        betaa[1] <- sample_l
        jieguo_old <- cbind(cpg, betaa)
        jieguo_old <- as.matrix(jieguo_old)
        for(i in 2:length(tcga_dir)) {
		#for(i in 121:130) {
		    message("file",file_num," is over")
		    file_num <- file_num+1
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

            file_l <- data.table::fread(file_name,header=FALSE)
            cpg <- file_l[,1]
            betaa <- file_l[,2]
            betaa[1] <- sample_l
			betaa <- as.vector(unlist(betaa))
            jieguo <- cbind(jieguo_old, betaa)
			jieguo_old <- jieguo
			rm(dirr_l,aa,file_name,sample_l,cpg,file_l,betaa,jieguo)
			gc()
			
        }
		
		
    file1 <- jieguo_old
    file2 <- file_ll
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


    #ventricle <- data.table::fread("methy_gene.txt",sep="\t",header=FALSE)
    #gene_ave(ventricle
    ave_result <- gene_ave(rep1_result)
	
    } else {message("please give your directory of methylation data!")}

}
