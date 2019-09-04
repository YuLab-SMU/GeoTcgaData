
#' Title
#'
#' @param file1 input file, a data.frame or a matrix
#' @param file2 output file name,a string
#' @param string a string,sep of the gene
#'
#' @return a matrix
#' @export
#'
#' @examples
#' aa <- c("MARCH1 /// MMA","MARC1","MARCH2 /// MARCH3","MARCH3 /// MARCH4","MARCH1")
#' bb <- c("2.969058399","4.722410064","8.165514853","8.24243893","8.60815086")
#' cc <- c("3.969058399","5.722410064","7.165514853","6.24243893","7.60815086")
#' file1 <- data.frame(aa=aa,bb=bb,cc=cc)

#' rep1(file1,"jieguo.txt"," /// ")

rep1<-function(file1,file2,string){
  file1<-as.matrix(file1)
  utils::write.table(t(as.matrix(colnames(file1))),file2,sep = "\t", append = T, row.names = F, 
                col.names = F, quote = F)
  for(i in 1 : dim(file1)[1]){
    gene<-unlist(strsplit(file1[i,1],string))
    for(j in 1 : length(gene)){
      utils::write.table(cbind(gene[j],matrix(file1[i,-1],nrow=1)),file2,sep="\t",append=T,
	  row.names=F,col.names=F,quote=F)
    }
}

}

#' Title
#'
#' @param file1 input file, a data.frame or a matrix
#' @param file2 output file name,a string
#' @param string a string,sep of the gene
#'
#' @return a matrix
#' @export
#'
#' @examples
#' aa <- c("MARCH1 /// MMA","MARC1","MARCH2 /// MARCH3","MARCH3 /// MARCH4","MARCH1")
#' bb <- c("2.969058399","4.722410064","8.165514853","8.24243893","8.60815086")
#' cc <- c("3.969058399","5.722410064","7.165514853","6.24243893","7.60815086")
#' file1 <- data.frame(aa=aa,bb=bb,cc=cc)
#' rep2(file1,"jieguo.txt"," /// ")
rep2<-function(file1,file2,string){
  file1<-as.matrix(file1)
  utils::write.table(t(as.matrix(colnames(file1))),file2,sep = "\t", append = T, row.names = F, 
                col.names = F, quote = F)
  a<-NULL
  for(i in 1 : dim(file1)[1]){
    if(length(grep(string,file1[i,1]))==0)
	{
	  a<-c(a,i)
	}
  }
  result<-file1[a,]
  utils::write.table(result,file2,sep="\t",append=T,row.names=F,col.names=F,quote=F)

}



