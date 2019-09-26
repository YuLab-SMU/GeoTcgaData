
#' Title preparer file for chi-square test
#'
#' @param jieguo2 result of ann_merge()
#'
#' @return a matrix
#' @export
#'
#' @examples
#' jieguo3 <- matrix(c(-1.09150,-1.47120,-0.87050,-0.50880,
#' -0.50880,2.0,2.0,2.0,2.0,2.0,2.601962,2.621332,2.621332,
#' 2.621332,2.621332,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,
#' 2.0,2.0,2.0,2.0,2.0,2.0,2.0),nrow=5)
#' rownames(jieguo3) <- c("AJAP1","FHAD1","CLCNKB","CROCCP2","AL137798.3")
#' colnames(jieguo3) <- c("TCGA-DD-A4NS-10A-01D-A30U-01","TCGA-ED-A82E-01A-11D-A34Y-01",
#' "TCGA-WQ-A9G7-01A-11D-A36W-01","TCGA-DD-AADN-01A-11D-A40Q-01",
#' "TCGA-ZS-A9CD-10A-01D-A36Z-01","TCGA-DD-A1EB-11A-11D-A12Y-01")
#' cnv_chi_file <- prepare_chi(jieguo3)
prepare_chi <- function(jieguo2) {

file1<-jieguo2
samples<-colnames(file1)
sampless<-samples
hehe<-samples
for(i in 1:length(samples))
{
	aa<-unlist(strsplit(samples[i],"-"))[4]
	bb<-unlist(strsplit(aa,""))[1]
	sampless[i]<-bb
	hehe[i]<-aa
}

cnv_chi <- matrix(c("gene","normalCNV","normalWild","tumorCNV","tumorWild"),nrow=1)
for(i in seq_len(dim(file1)[1]))
{
	normalCNV<-0
	normalWild<-0
	tumorCNV<-0
	tumorWild<-0
	for(j in 2:dim(file1)[2])
	{
		if((sampless[j-1]=="1")&&(abs(as.numeric(file1[i,j])-2)>0.5))
		{
			normalCNV<-normalCNV+1
		}
		if((sampless[j-1]=="1")&&(abs(as.numeric(file1[i,j])-2)<=0.5))
		{
			normalWild<-normalWild+1
		}
		if((sampless[j-1]=="0")&&(abs(as.numeric(file1[i,j])-2)>0.5))
		{
			tumorCNV<-tumorCNV+1
		}
		if((sampless[j-1]=="0")&&(abs(as.numeric(file1[i,j])-2)<=0.5))
		{
			tumorWild<-tumorWild+1
		}
	}
	jieguo<-c(rownames(file1)[i],normalCNV,normalWild,tumorCNV,tumorWild)
	jieguo<-matrix(jieguo,nrow=1)
	cnv_chi <- rbind(cnv_chi,jieguo)

}
return(cnv_chi)
}


#' Title do chi-square test to find differential genes
#'
#' @param rt result of prepare_chi()
#'
#' @return a matrix
#' @export
#'
#' @examples
#' jieguo3 <- matrix(c(-1.09150,-1.47120,-0.87050,-0.50880,
#'                     -0.50880,2.0,2.0,2.0,2.0,2.0,2.601962,2.621332,2.621332,
#'                     2.621332,2.621332,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,
#'                     2.0,2.0,2.0,2.0,2.0,2.0,2.0),nrow=5)
#' rownames(jieguo3) <- c("AJAP1","FHAD1","CLCNKB","CROCCP2","AL137798.3")
#' colnames(jieguo3) <- c("TCGA-DD-A4NS-10A-01D-A30U-01","TCGA-ED-A82E-01A-11D-A34Y-01",
#' "TCGA-WQ-A9G7-01A-11D-A36W-01","TCGA-DD-AADN-01A-11D-A40Q-01",
#' "TCGA-ZS-A9CD-10A-01D-A36Z-01","TCGA-DD-A1EB-11A-11D-A12Y-01")
#' rt <- prepare_chi(jieguo3)
#' chiResult <- differential_cnv(rt)
differential_cnv <- function(rt) {
    colnames(rt) <- rt[1,]
    rt <- rt[-1,]
    rownames(rt) <- rt[,1]
    rt <- rt[,-1]
    rtt <- matrix(as.numeric(rt),ncol=ncol(rt))
    rownames(rtt) <- rownames(rtt)
    outTab=data.frame()
    for(i in seq_len(dim(rtt)[1])){
        x=matrix(c(rtt[i,1],rtt[i,2],rtt[i,3],rtt[i,4]), ncol = 2)
        chiTest=stats::chisq.test(x)
        normalRatio=rtt[i,1]/(rtt[i,1]+rtt[i,2])
        tumorRatio=rtt[i,3]/(rtt[i,3]+rtt[i,4])
        Gene=row.names(rtt[i,])
        Stat=chiTest$statistic
        Pvalue=chiTest$p.value
        outTab=rbind(outTab,cbind(Gene,normalRatio,tumorRatio,Stat,Pvalue))
    }
    pvalue=as.numeric(as.vector(outTab[,"Pvalue"]))
    adjP=stats::p.adjust(pvalue,method ="bonferroni")
    outTab=cbind(outTab,adjPvalue=adjP)
	rownames(outTab) <- rownames(rt)
	return(outTab)
}




