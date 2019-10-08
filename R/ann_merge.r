
deal_meta <- function(metadata_file) {
    file1<-as.matrix(utils::read.table(metadata_file,sep="\t",header=F))

    file2<-file1[grep(pattern="entity_id",file1[,1]),]
    file3<-file1[grep(pattern="entity_submitter_id",file1[,1]),]

    file3<-file3[nchar(file3)==57]

    file2<-gsub("      entity_id: ","",file2)
    file2<-gsub(", ","",file2)
    file3<-gsub("      entity_submitter_id: ","",file3)
    file3<-gsub(", ","",file3)
	return(cbind(file2,file3))
}

#' Merge the copy number variation data downloaded from TCGA using gdc
#'
#' @param dirr a string of direction, catalogue of copy number variation data
#' @param metadatafile a metadata file download from TCGA
#'
#' @return a matrix,each column is a sample, each row is a gene
#' @export
#'
#' @examples
#' metadatafile_name <- "metadata.cart.2018-11-09.json"
#' \dontrun{jieguo2 <- ann_merge(dirr = system.file(file.path("extdata","cnv"),
#' package="GeoTcgaData"),metadatafile=metadatafile_name)}
ann_merge <- function(dirr,metadatafile) {
    tcga_dir <- dir(dirr)
	temp_dirr <- tempdir()
	if(!file.exists(file.path(temp_dirr,"temp"))) {
	  dir.create(file.path(temp_dirr,"temp"))
	}
	temp_dir <- file.path(temp_dirr,"temp")
	lf<-list.files(temp_dir)
    file.remove(file.path(temp_dir,lf))
	for(k in seq_len(length(tcga_dir))) {
    dirr_l <- file.path(dirr, tcga_dir[k])
    files <- dir(dirr_l)


    for(j in seq_len(length(files))) {
        if(length(grep(".txt",files[j]))>0) {

            file_name <- file.path(dirr_l, dir(dirr_l)[j])
        }
    }

    aa <- as.matrix(data.table::fread(file_name,header=TRUE))
    aa<-aa[which(abs(as.numeric(aa[,6]))>0.2),]
    genes<-rep("he",dim(aa)[1])
	    for(i in seq_len(dim(aa)[1])) {

            for(j in seq_len(dim(genePos)[1])) {
                if(aa[i,2]==genePos[j,5]) {
        	        if((as.numeric(aa[i,3])>as.numeric(genePos[j,3]))|| (as.numeric(aa[i,4])<as.numeric(genePos[j,2]))){}
        	        else {
        		        genes[i]<-paste(genes[i],genePos[j,1],sep=",")
        	       }
                }

            }

	    }

	bb<-cbind(aa,genes)
	bb<-bb[which(bb[,7]!="he"),]
#when bb is just one line, do this
    if(length(bb)<10)
    {
        bb[7]<-gsub("he,","",bb[7])


        dd<-unlist(strsplit(bb[7],","))
        for(j in 1:length(dd))
          {
            jieguo<-c(bb[1:6],dd[j])
            jieguo<-matrix(jieguo,nrow=1)
            utils::write.table(jieguo,file.path(temp_dir,paste(aa[1,1],".dada",sep="")),sep="\t",row.names=F,col.names=F,append=T,quote=F)
         }


#when bb is more than one line, do this
    }else{
        bb[,7]<-gsub("he,","",bb[,7])
            for(i in 1:dim(bb)[1])
        {
        dd<-unlist(strsplit(bb[i,7],","))
        for(j in 1:length(dd))
          {
            jieguo<-c(bb[i,1:6],dd[j])
            jieguo<-matrix(jieguo,nrow=1)
            utils::write.table(jieguo,file.path(temp_dir,paste(aa[1,1],".dada",sep="")),sep="\t",row.names=F,col.names=F,append=T,quote=F)
         }
       }

    }


}

filenames<-dir(temp_dir)
genes<-NULL
file_id<-filenames
i=1
for(file in filenames)
{
	aa<-as.matrix(utils::read.table(file.path(temp_dir,file),sep="\t",header=F))

	genes<-union(genes,aa[,7])
	file_id[i]<-aa[1,1]
	i=i+1

}

jieguo<-matrix(0,length(genes),length(filenames))
rownames(jieguo)<-genes
colnames(jieguo)<-file_id

for(file in filenames)
{
	aa<-as.matrix(utils::read.table(file.path(temp_dir,file),sep="\t",header=F))
	for(i in 1:dim(aa)[1])
	{
		#jieguo[aa[i,7],aa[1,1]]=2^(as.numeric(aa[i,6])+1)
		jieguo[aa[i,7],aa[1,1]]=aa[i,6]

	}
}

metadata<-deal_meta(metadatafile)
rownames(metadata)<-metadata[,1]
filenames_tcga<-metadata[file_id,2]
colnames(jieguo)<-filenames_tcga
#jieguo<-cbind(genes,jieguo)
#jieguo<-rbind(c("genes",filenames_tcga),jieguo)

for(i in 1:dim(jieguo)[1])
{
	for(j in 2:dim(jieguo)[2])
	{
		jieguo[i,j]<-2^(1+as.numeric(jieguo[i,j]))
	}
}
lf<-list.files(temp_dir)
file.remove(file.path(temp_dir,lf))

return(jieguo)

}


