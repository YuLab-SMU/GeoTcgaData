
#' Find the mean value of the gene in each module
#'
#' @param geneExpress a data.frame
#' @param module a data.frame
#'
#' @return a matrix, means the mean of gene expression value in
#' the same module
#' @export
#'
#' @examples
#' result <- cal_mean_module(geneExpress,module)
cal_mean_module<-function(geneExpress,module){
  out <- file.path(tempdir(), "module_result.txt")
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
	utils::write.table(cbind(module[i,1],module_mean),out,sep="\t",
	append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
  output_module <- data.table::fread(out,sep="\t",header=FALSE)
  file.remove(out)
  return(output_module)

}

