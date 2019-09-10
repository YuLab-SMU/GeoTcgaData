

#' Convert  ENSEMBL gene id to gene Symbol in TCGA
#'
#' @param profile a data.frame
#'
#' @return a data.frame, gene symbols and their expression value
#' @export
#'
#' @examples
#' result <- id_conversion(profile)
id_conversion<-function(profile){
  file2<-as.matrix(profile)
  for(i in 1:dim(profile)[1]){
    file2[i,1]<-unlist(strsplit(file2[i,1],"\\."))[1]
  }


  file3<-hgnc
  rownames(file3)<-file3[,5]
  haha<-intersect(file2[,1],file3[,5])
  hehe<-file3[haha,c(1,5)]
  rownames(file2)<-file2[,1]
  file4<-file2[hehe[,2],]
  file4[,1]<-hehe[,1]
  return(file4)
}

