
#' Title
#'
#' @param geneExpress a data.frame
#' @param module a data.frame
#' @param result a string
#'
#' @return a matrix
#' @export
#'
#' @examples
#' cal_mean_module(geneExpress,module,"result.txt")
cal_mean_module<-function(geneExpress,module,result){
  geneExpress=as.matrix(geneExpress)
  rownames(geneExpress)<-geneExpress[,1]
  genes<-geneExpress[,1]
  geneExpress<-geneExpress[,-1]
  for(i in 1:dim(module)[1]){
    modulen<-unlist(strsplit(module[i,2],","))
    modulen<-intersect(modulen,genes)
    n<-length(modulen)
    modulen_matrix<-geneExpress[modulen[1],]
    for(j in 2:n){

      modulen_matrix<-rbind(modulen_matrix,geneExpress[modulen[j],])
    }
    modulen_matrix<-matrix(as.numeric(modulen_matrix),nrow<-nrow(modulen_matrix))
    module_mean<-colMeans(modulen_matrix)
    module_mean<-matrix(module_mean,nrow<-1)
	utils::write.table(cbind(module[i,1],module_mean),result,sep="\t",
	append=T,row.names=F,col.names=F,quote=F)
  }


}

