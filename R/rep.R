
#' Handle the case where one id corresponds to multiple genes
#'
#' @param input_file1 input file, a data.frame or a matrix
#' @param string a string,sep of the gene
#'
#' @return a data.frame, when an id corresponds to multiple genes, 
#' the expression value is assigned to each gene
#' @export
#'
#' @examples
#' aa <- c("MARCH1 /// MMA","MARC1","MARCH2 /// MARCH3","MARCH3 /// MARCH4","MARCH1")
#' bb <- c("2.969058399","4.722410064","8.165514853","8.24243893","8.60815086")
#' cc <- c("3.969058399","5.722410064","7.165514853","6.24243893","7.60815086")
#' input_fil <- data.frame(aa=aa,bb=bb,cc=cc)

#' rep1_result <- rep1(input_fil," /// ")

rep1<-function(input_file1,string){
  input_file1 = as.matrix(input_file1)
  out <- file.path(tempdir(), "rep1_result.txt")
  utils::write.table(t(as.matrix(colnames(input_file1))),out,sep = "\t", append = TRUE, row.names = FALSE, 
                col.names = FALSE, quote = FALSE)
  for(i in 1 : dim(input_file1)[1]){
    gene<-unlist(strsplit(input_file1[i,1],string))
    for(j in 1 : length(gene)){
      utils::write.table(cbind(gene[j],matrix(input_file1[i,-1],nrow=1)),out,sep="\t",append=TRUE,
	  row.names=FALSE,col.names=FALSE,quote=FALSE)
    }
}
  output_rep1 <- data.table::fread(out,sep="\t",header=F)
  file.remove(out)
  return(output_rep1)

}

#' Handle the case where one id corresponds to multiple genes
#'
#' @param input_file1 input file, a data.frame or a matrix
#' @param string a string,sep of the gene
#'
#' @return a matrix,when an id corresponds to multiple genes,
#' the expression value is deleted
#' @export
#'
#' @examples
#' aa <- c("MARCH1 /// MMA","MARC1","MARCH2 /// MARCH3","MARCH3 /// MARCH4","MARCH1")
#' bb <- c("2.969058399","4.722410064","8.165514853","8.24243893","8.60815086")
#' cc <- c("3.969058399","5.722410064","7.165514853","6.24243893","7.60815086")
#' input_fil <- data.frame(aa=aa,bb=bb,cc=cc)
#' rep2_result <- rep2(input_fil," /// ")
rep2<-function(input_file1,string){
  input_file1 = as.matrix(input_file1)

  a<-NULL
  for(i in 1 : dim(input_file1)[1]){
    if(length(grep(string,input_file1[i,1]))==0)
	{
	  a<-c(a,i)
	}
  }
  result<-input_file1[a,]
  return(result)
}



