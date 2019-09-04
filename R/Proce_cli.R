
#' Title
#'
#' @param Files_dir1 a dir data
#'
#' @return a matrix
#' @export
#'
#' @examples
#' tcga_cli_deal("your_clinical_directory")
tcga_cli_deal<-function(Files_dir1){
  if(Files_dir1!="your_clinical_directory") {
  Files_dir <- dir(Files_dir1)
  Files=Files_dir
  file_id<-Files
  for(i in 1:length(file_id)){
    file_id[i]<-unlist(strsplit(file_id[i],"\\."))[3]
  }
  time<-"</clin_shared:days_to_last_followup>"
  state<-"</clin_shared:vital_status>"

  days="</clin_shared:days_to_death>"
  for(file in Files){
    #haha means the sample has no survival data
    timee="haha"
    statee<-"haha"
    aa<-as.matrix(utils::read.table(file,sep="\t",header=T))
    for(i in 1:dim(aa)[1]){
        if(length(grep(time,aa[i,1]))>0){
             timee<-unlist(strsplit(aa[i,1],">"))[2]

             timee<-gsub("</clin_shared:days_to_last_followup","",timee)
        }
    }

    for(i in 1:dim(aa)[1]){
        if(length(grep(state,aa[i,1]))>0){
             statee<-unlist(strsplit(aa[i,1],">"))[2]

             statee<-gsub("</clin_shared:vital_status","",statee)
        }
    }
        for(i in 1:dim(aa)[1]){
        if(length(grep(days,aa[i,1]))>0){
             timee<-unlist(strsplit(aa[i,1],">"))[2]

             timee<-gsub("</clin_shared:days_to_death","",timee)
        }
    }
    file_idd<-unlist(strsplit(file,"\\."))[3]

    utils::write.table(cbind(file_idd,timee,statee),"clin2.txt",
	sep="\t",row.names=F,col.names=F,append=T,quote=F)
  }
} else {message("please give your directory!")} }
