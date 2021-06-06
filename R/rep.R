#' Handle the case where one id corresponds to multiple genes
#'
#' @param input_file input file, a data.frame or a matrix
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
#' input_file <- data.frame(aa=aa,bb=bb,cc=cc)

#' rep1_result <- rep1(input_file," /// ")

rep1 <- function(input_file,string){
    name <- colnames(input_file)[1]
    genelist <- strsplit(input_file[,1], string)
    geneLength <- unlist(lapply(genelist, length))
    input_file <- input_file[, -1]
    output <- apply(input_file, 2, rep, times = geneLength)
    output2 <- matrix(as.numeric(output), nrow = nrow(output))
    colnames(output2) <- colnames(output)
    output2 <- data.frame(unlist(genelist), output2, check.names = FALSE)
    colnames(output2)[1] <- name
    output2
}

#' Handle the case where one id corresponds to multiple genes
#'
#' @param input_file input file, a data.frame or a matrix
#' @param string a string,sep of the gene
#'
#' @return a data.frame, when an id corresponds to multiple genes,
#' the expression value is deleted
#' @export
#'
#' @examples
#' aa <- c("MARCH1 /// MMA","MARC1","MARCH2 /// MARCH3","MARCH3 /// MARCH4","MARCH1")
#' bb <- c("2.969058399","4.722410064","8.165514853","8.24243893","8.60815086")
#' cc <- c("3.969058399","5.722410064","7.165514853","6.24243893","7.60815086")
#' input_file <- data.frame(aa=aa,bb=bb,cc=cc)
#' rep2_result <- rep2(input_file," /// ")
rep2 <- function(input_file, string){
  unKeep <- grep(string, input_file[,1])
  if (length(unKeep) > 0) input_file <- input_file[-unKeep, ]
  input_file
}



