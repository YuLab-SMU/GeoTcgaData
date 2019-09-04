
#' Title
#'
#' @param file1 a data.frame
#' @param k a number
#' @param result a string for output file name
#'
#' @return a data.frame
#' @export
#'
#' @examples
#' aa <- c("Gene Symbol","MARCH1","MARC1","MARCH1","MARCH1","MARCH1")
#' bb <- c("GSM1629982","2.969058399","4.722410064","8.165514853","8.24243893","8.60815086")
#' cc <- c("GSM1629982","3.969058399","5.722410064","7.165514853","6.24243893","7.60815086")
#' file1 <- data.frame(aa=aa,bb=bb,cc=cc)
#' gene_ave(file1)
gene_ave<-function(file1,k=1,result="result.txt"){
  file1<-as.matrix(file1)
  rownames(file1)<-file1[,k]
  x<-file1
  ID<-rownames(file1)
  ID <- factor(ID,levels=unique(ID))
#	rowsum(x,ID,reorder=FALSE,na.rm=TRUE)/as.vector(table(ID))
    
  x<-as.data.frame(x)
  for(gg in 2:dim(x)[2]){
	x[,gg]<-as.numeric(as.character(x[,gg]))
  }
  aa<-x[,-k]
  y <- rowsum(aa,ID,reorder=FALSE,na.rm=TRUE)
  n <- rowsum(1L-is.na(aa),ID,reorder=FALSE)
  gege<-y/n
  gege[1,] <- file1[1,-1]
  utils::write.table(gege,result,sep="\t",row.names=T,col.names=F,append=T,quote=F)
} 











